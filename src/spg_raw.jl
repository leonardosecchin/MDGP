###############################
# SPG for Euclidean coordinates
# (raw stress)
###############################
include("spg_core.jl")

# point/vector structure
struct SPG_VECTOR_EUC_RAW
    X::Matrix{Float64}
    d::Vector{Float64}
end


# Operations between SPG_VECTOR_EUC_RAWs

# Z = X
function spg_vector_cp!(
    Z::SPG_VECTOR_EUC_RAW, X::SPG_VECTOR_EUC_RAW, work::SPG_WORKSPACE
)
    Z.X .= X.X
    Z.d .= X.d
end

# dot(X,Y)
function spg_vector_dot(
    X::SPG_VECTOR_EUC_RAW, Y::SPG_VECTOR_EUC_RAW, work::SPG_WORKSPACE
)
    r = dot(X.X, Y.X)
    @inbounds for k in work.idxD
        r += X.d[k] * Y.d[k]
    end
    return r
end

# Z = X + a*Y
function spg_vector_xpay!(
    Z::SPG_VECTOR_EUC_RAW,
    X::SPG_VECTOR_EUC_RAW,
    a::Float64,
    Y::SPG_VECTOR_EUC_RAW,
    work::SPG_WORKSPACE
)
    @. Z.X = X.X + a * Y.X
    @inbounds for k in work.idxD
        Z.d[k] = X.d[k] + a * Y.d[k]
    end
end

# supnorm(X)
function spg_vector_supn(X::SPG_VECTOR_EUC_RAW, work::SPG_WORKSPACE)
    r = norm(X.X, Inf)
    @inbounds for k in work.idxD
        r = max(r, abs(X.d[k]))
    end
    return r
end

# projected direction
function proj_d!(
    data::DATA, d::SPG_VECTOR_EUC_RAW, lambda::Float64, work::SPG_WORKSPACE
)
    @. d.X = -lambda * work.g.X
    @inbounds for k in work.idxD
        d.d[k] = -lambda * work.g.d[k]
    end

    @inbounds @views for k in work.idxD
        L, U = data.D[k,1], data.D[k,2]
        xk = work.x.d[k]
        aux = xk + d.d[k]
        if aux < L
            d.d[k] = L - xk
        elseif aux > U
            d.d[k] = U - xk
        end
    end
end


# Stress function
#
# stress(x) = 0.5 * sum_{k-th distance} w[k] * (||x.X_i - x.X_j|| - x.d[k])^2
#
# where i,j are the vertices corresponding to the k-th distance
function stress(
    data::DATA,
    x::SPG_VECTOR_EUC_RAW,
    work::SPG_WORKSPACE
)
    σ = 0.0

    @inbounds @views for k in work.idxD
        i, j = data.Dij[k,1], data.Dij[k,2]
        work.dists[k] = d(i, j, x.X)
        σ += work.w[k] * (work.dists[k] - x.d[k])^2
    end
    if isnan(σ)
        error("stress returned NaN")
    end

    return σ/2.0
end


# gradient of stress
function grad_stress!(
    data::DATA,
    x::SPG_VECTOR_EUC_RAW,
    g::SPG_VECTOR_EUC_RAW,
    work::SPG_WORKSPACE
)
    work.work .= 0.0
    g.X .= 0.0

    @inbounds @views for k in work.idxD
        wk = work.w[k]
        dk = x.d[k]
        distk = work.dists[k]

        if distk > 0.0
            i, j = data.Dij[k,1], data.Dij[k,2]
            aux = wk * (1.0 - dk/distk)
            work.work[i] += aux
            work.work[j] += aux
            @. g.X[1:3,i] -= aux * x.X[1:3,j]
            @. g.X[1:3,j] -= aux * x.X[1:3,i]
        end

        # gradient w.r.t. d
        g.d[k] = wk * (dk - distk)
    end

    @inbounds @views for i in 1:data.nv
        @. g.X[1:3,i] += work.work[i] * x.X[1:3,i]
    end

    return
end
