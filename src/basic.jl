struct WORKSPACE
    U::Matrix{Float64}
    w::Vector{Float64}
    dists::Vector{Float64}
    Zd::Vector{Float64}
    GX::Matrix{Float64}
    Gd::Vector{Float64}
    Pd::Vector{Float64}
    ZXnew::Matrix{Float64}
    Zdnew::Vector{Float64}
    GXnew::Matrix{Float64}
    Gdnew::Vector{Float64}
    SX::Matrix{Float64}
    YX::Matrix{Float64}
    Sd::Vector{Float64}
    Yd::Vector{Float64}
    ZXbest::Matrix{Float64}
    Zdbest::Vector{Float64}
    idxD::Matrix{Int64}
    work::Vector{Float64}
    lastf::Vector{Float64}
end

function init_workspace(nv, nd, lsm)
    return WORKSPACE(
        Matrix{Float64}(undef,3,3),  # U
        Vector{Float64}(undef,nd),   # w
        Vector{Float64}(undef,nd),   # dists
        Vector{Float64}(undef,nd),   # Zd
        Matrix{Float64}(undef,3,nv), # GX
        Vector{Float64}(undef,nd),   # Gd
        Vector{Float64}(undef,nd),   # Pd
        Matrix{Float64}(undef,3,nv), # ZXnew
        Vector{Float64}(undef,nd),   # Zdnew
        Matrix{Float64}(undef,3,nv), # GXnew
        Vector{Float64}(undef,nd),   # Gdnew
        Matrix{Float64}(undef,3,nv), # SX
        Matrix{Float64}(undef,3,nv), # YX
        Vector{Float64}(undef,nd),   # Sd
        Vector{Float64}(undef,nd),   # Yd
        Matrix{Float64}(undef,3,nv), # ZXbest
        Vector{Float64}(undef,nd),   # Zdbest
        Matrix{Int64}(undef,nv,2),   # idxD
        Vector{Float64}(undef,nv),   # work
        fill(-Inf, lsm)              # lastf
    )
end

# return [ lower d(i,j) ; upper d(i,j)] from D. Distance d(i,j) must exists in D
@inline function d(i, j, D, ij_to_D)
    return @inbounds @views D[ij_to_D[i,j], 1:end]
end


# bond length  : r_i = d_{i-1,i}
# bond angle   : theta_i = "angle between segments i-2 --- i-1 and i-1 --- i"
# torsion angle: omega_i = "angle between the normals through the planes deﬁned by the atoms i-3, i-2, i-1 and i-2, i-1, i"

# given vertices i-2, i-1, i, compute the cosine of the bond angle theta_i
# by the law of cosines:
#
# cos(theta_i) = d_{i-1,i}^2 + d_{i-2,i-1}^2 - d_{i-2,i}^2
#                -----------------------------------------
#                         2 d_{i-1,i} d_{i-2,i-1}

function costheta(d10, d20, d21)
    cosine = (d10^2 + d21^2 - d20^2)/(2.0*d21*d10)
    return sign(cosine) * min(abs(cosine),1.0)
end


# Given vertices i-3, i-2, i-1, i, compute the cosine of the torsion angle omega_i:
#
# cos(omega_i) = r^2_{i-2} + d^2_{i-2,i} - 2 r_{i-2} d_{i-2,i} cos(theta_{i-1}) cos(theta_i) - d^2_{i-3,i}
#                -----------------------------------------------------------------------------------------
#                                    2 r_{i-2} d_{i-2,i} sin(theta_{i-1}) sin(theta_i)
#
#              = costheta(d32,d30,d20) - cos(theta_{i-1}) cos(theta_i)
#                -----------------------------------------------------
#                            sin(theta_{i-1}) sin(theta_i)
#
# The (possible) interval distance d30 = d^2_{i-3,i} must be fixed outside this function.
# see Lavor, Liberti, Maculan. A note on "A branch-and-prune algorithm for the molecular distance geometry problem". Intl. Trans. in Op. Res. 18 (2011) 751-752

function cosomega(d10, d32, d31, d21, d20, d30)
    costheta1 = costheta(d21,d31,d32)
    costheta0 = costheta(d20,d10,d21)
    r         = costheta(d32,d30,d20)

    cosine = (r - costheta1 * costheta0) / (sqrt(1.0 - costheta1^2) * sqrt(1.0 - costheta0^2))
    return sign(cosine) * min(abs(cosine),1.0)
end


# Compute the distance d_{i-3,i} given all exact distances between vertices i-3, i-2, i-1, i and the cosine of the torsion angle omega_i

function d30(d10, d32, d31, d21, d20, cosomega)
    costheta1 = costheta(d32,d31,d21)
    costheta0 = costheta(d20,d10,d21)

    sintheta0 = sqrt(1.0 - costheta0^2)
    sintheta1 = sqrt(1.0 - costheta1^2)

    # isolate d_{i-3,i} in the formula for the cosine of omega
    return sqrt(d32^2 + d20^2 - 2.0*d32*d20*(cosomega*sintheta0*sintheta1 + costheta0*costheta1))
end


# Mean/Largest Distance Error (MDE/LDE)
function LDE(Dij, D, idxD, X::Matrix{Float64})
    lde = 0.0

    @inbounds @views for k in idxD
        dist = euclidean(X[1:3,Dij[k,1]], X[1:3,Dij[k,2]])
        L,U = D[k,1:2]
        if L == U
            lde = max(lde, abs(L - dist)/L)
        else
            lde = max(lde, max( max(L - dist, 0.0)/L, max(dist - U, 0.0)/U ) )
        end
    end

    return lde
end

function MDE_LDE(Dij, D, idxD, X::Matrix{Float64})
    mde = 0.0
    lde = 0.0

    @inbounds @views for k in idxD
        dist = euclidean(X[1:3,Dij[k,1]], X[1:3,Dij[k,2]])
        L,U = D[k,1:2]
        if L == U
            mde += abs(L - dist)/L
            lde = max(lde, abs(L - dist)/L)
        else
            mde += max(L - dist, 0.0)/L + max(dist - U, 0.0)/U
            lde = max(lde, max( max(L - dist, 0.0)/L, max(dist - U, 0.0)/U ) )
        end
    end

    return mde/length(idxD), lde
end

MDE_LDE(Dij, D, X::Matrix{Float64}) = MDE_LDE(Dij, D, 1:size(D,1), X)

function MDE_LDE(Dij, D, idxD, dists::Vector{Float64})
    mde = 0.0
    lde = 0.0

    @inbounds @views for k in idxD
        L,U = D[k,1:2]
        if L == U
            mde += abs(L - dists[k])/L
            lde = max(lde, abs(L - dists[k])/L)
        else
            mde += max(L - dists[k], 0.0)/L + max(dists[k] - U, 0.0)/U
            lde = max(lde, max( max(L - dists[k], 0.0)/L, max(dists[k] - U, 0.0)/U ) )
        end
    end

    return mde/length(idxD), lde
end

# centralize conformation w.r.t. all atoms
function centralize!(X, Xout, idx)
    @inbounds @views work = sum(X[1:3,idx], dims=2)./length(idx)
    @inbounds @views for i in idx
        @. Xout[1:3,i] = X[1:3,i] - work
    end
end


# RMSD between two conformations w.r.t. given atoms indices,
# conformations are considered centralized
function rmsd_protein(X1, X2, idx)
    @inbounds @views SVD = svd(X2[1:3,idx] * X1[1:3,idx]')
    return @inbounds @views (1.0/sqrt(length(idx))) * norm(SVD.U * SVD.Vt * X1[1:3,idx] - X2[1:3,idx], 2)
end
