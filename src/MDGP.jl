"""
Multistart strategy for Molecular Distance Geometry Problem
"""

module MDGP

using LinearAlgebra
using SparseArrays
using Printf
using Distances
using Random
using Distributions
using DelimitedFiles

export MDGP_multistart, MDGP_read

include("vectors.jl")
include("basic.jl")
include("preprocess.jl")
include("spg.jl")
include("conformation.jl")
include("read.jl")

########################
# Main function
########################
"""
    sols, types, ldes, mdes, time_init, time_total = MDGP_multistart(Dij, D, P, atoms, res, torsions; [OPTIONS])

Multistart strategy for MDGP.
"""
function MDGP_multistart(
            Dij_orig::Matrix{Int64},         # nd x 2 matrix of indices of distances
            D_orig::Matrix{Float64},         # nd x 2 matrix of distances
            P_orig::Matrix{Int64},           # predecessors and branching signs
            atoms::Vector{String},           # atoms names
            res::Vector{Int64},              # residues indices
            torsions::Matrix{Float64};       # torsion angles and left/right displacements (in degrees)
            # SPG options
            spg_maxit::Int64     = 30000,    # maximum number of SPG iterations
            spg_lacktol::Real    = 1e-8,     # tolerance to declare lack of progress in SPG
            spg_eta::Real        = 1e-4,     # Armijo's parameter
            spg_lsm::Int64       = 10,       # length of the history for non-monotone line search
            spg_lmin::Real       = 1e-20,    # minimum value for spectral steplength
            spg_lmax::Real       = 1e+20,    # maximum value for spectral steplength
            # general options
            num_X0::Int64        = 50,       # number of initial conformations
            num_impr::Int64      = 5,        # number of improvement trials
            max_init::Int64      = 500,      # max number of initial trials
            max_torsion::Int64   = 20,       # max number of torsion angles trials
            max_sols::Int64      = 1,        # number of solutions required
            max_time::Real       = 7200,     # max time in seconds
            max_similar::Int64   = 50,       # max number of consecutive similar init conf
            # tolerances
            tol_lde::Real        = 1e-3,     # tolerance for optimality (LDE)
            tol_mde::Real        = 1e-3,     # tolerance for optimality (MDE)
            tol_stress::Real     = 1e-7,     # tolerance for optimality (SPG, stress)
            tol_exact::Real      = 1e-12,    # maximum interval length to consider a distance exact
            tol_similar::Real    = 5.0,      # tolerance to consider two conformations equal (RMSD)
            # other
            seed::Int64          = -1,       # random seed (<0 for any)
            verbose::Int64       = 1         # output level
            )

    @assert num_X0 > 0 "num_X0 must be positive"
    @assert max_init > 0 "max_init must be positive"
    @assert max_torsion >= 0 "max_torsion must be non negative"
    @assert max_sols > 0 "max_sols must be positive"
    @assert max_time > 0.0 "max_time must be positive"

    @assert spg_lsm > 0 "spg_lsm must be positive"
    @assert (spg_eta > 0) && (spg_eta < 1.0) "spg_eta must be in (0,1)"
    @assert spg_lmin > 0.0 "spg_lmin must be positive"
    @assert spg_lmax > 0.0 "spg_lmax must be positive"

    @assert (tol_mde > 0.0) || (tol_lde > 0.0) "tol_lde or tol_lde must be positive"
    @assert tol_exact >= 0.0 "tol_exact must be positive"
    @assert tol_similar >= 0.0 "tol_similar must be non negative"

    if seed >= 0
        Random.seed!(seed)
    end

    status = -1

    time_pre = @elapsed begin

    Dij = deepcopy(Dij_orig)
    D = deepcopy(D_orig)

    # P has 4 columns: predecessors (cols 1 to 3) and branching signs (col 4)
    P = init_P(P_orig, maximum(Dij))

    if check_basics(Dij, D, P) < 0
        return [], [], [], [], Inf, Inf
    end

    # Consolidate repeated distances
    consolidate_distances!(Dij, D)

    # number of vertices and distances
    nv = maximum(Dij)
    nd = size(Dij,1)

    # construct the map (i,j) to the position in D
    ij_to_D = spzeros(Int64, nv, nv)
    @inbounds @views for k in 1:nd
        ij_to_D[ Dij[k,1], Dij[k,2] ] = k
    end
    ij_to_D = Symmetric(ij_to_D, :U)

    # check the existence of necessary distances
    if !check_necessarydistances(nv, D, P, ij_to_D, tol_exact)
        return [], [], [], [], Inf, Inf
    end

    # Try to improve d_{i-3,i} by computing it theoretically
    flag = tight_bounds!(nv, D, P, ij_to_D, tol_exact, verbose)
    if flag == 1
        @error "Problem is infeasible!"
        return [], [], [], [], Inf, Inf
    elseif flag == 2
        @warn "Instance is exact, enumerative strategies are recommended."
        return [], [], [], [], Inf, Inf
    end

    # groups of distances
    idxDpred = Int64[]      # between predecessors
    idxDnonpred = Int64[]   # extra
    idxDvdw = Int64[]       # Van der Walls
    @inbounds for k in 1:size(D,1)
        i,j = Dij[k,1:2]
        if (j in P[i,1:3]) || (i in P[j,1:3])
            push!(idxDpred, k)
        elseif D[k,2] < 900.0
            push!(idxDnonpred, k)
        else
            push!(idxDvdw, k)
        end
    end

    idxDpred = consec_range(idxDpred)
    idxDnonpred = consec_range(idxDnonpred)
    idxDvdw = consec_range(idxDvdw)

    idxDprednonpred = consec_range(union(idxDpred,idxDnonpred))

    # adjacent distances to each atom
    adj = Vector{Vector{Int64}}(undef, nv)
    for v in 1:nv
        adj[v] = Int64[]
        for k in 1:nd
            if Dij[k,2] == v
                push!(adj[v], k)
            end
        end
    end

    # ======================
    # First definitions and allocation
    # ======================

    X = start_conformation(nv, D, ij_to_D)

    # initialize workspace
    work = init_workspace(nv, nd, spg_lsm)

    # stress weights
    work.w .= 1.0
    @views work.w[idxDpred] .= 2.0      # discretization distances
    work.w ./= norm(work.w)

    stop = false

    total_count = 0

    fixed_torsions = Vector{Float64}(undef, nv)
    fixed_torsions[1:3] .= 0.0

    # vector of solutions
    sols = Vector{Matrix{Float64}}(undef, 0)
    ldes = Float64[]
    mdes = Float64[]
    soltypes = Int64[]
    stress = Float64[]

    strs = 0.0

    # similarity between generated conformations
    consec_similar = 0
    if nv <= 500
        rmsd_idx = 1:nv
    else
        rmsd_idx = (1:nv)[atoms .== "CA"]
    end
    Xtmp1 = similar(X)
    Xtmp2 = similar(X)

    end    # end of @elapsed begin

    if verbose > 0
        ndexact = count(D[:,2] .- D[:,1] .<= tol_exact)
        @printf("Time preprocess phase and allocation: %.6lf s\n", time_pre)
        println("Preprocessed data has $(nv) atoms, $(ndexact) exact distances and $(nd - ndexact) interval distances.")
    end

    # ======================
    # Starting conformations
    # ======================

    if verbose > 0
        println("\nComputing starting conformations...")
    end

    time_total = 0.0
    time_init = 0.0

    while (length(sols) < num_X0) && (total_count < max_init)

        time_initk = @elapsed begin

        if verbose > 0
            print("\rAttempt $(total_count), $(length(sols)) conformations found")
        end

        if count(soltypes .>= 0) >= max_sols
            stop = true
            break
        end

        # sort torsion angles
        @inbounds for v in 4:nv
            sgn = (P[v,4] == 0) ? rand([-1,1]) : P[v,4]
            fixed_torsions[v] = sort_torsion_angle(v, P, torsions, sgn)
        end

        # new conformation
        @inbounds if construct_conformation!(4,
                                            Dij, D, P, ij_to_D, torsions,
                                            X, fixed_torsions, adj, work,
                                            max_torsion, tol_lde, false
                                            )

            total_count += 1
            # Try to improve the new conformation by minimizing LDE
            # Note that we can neglect discretization distances as
            # they will still be satisfied.
            if num_impr > 0
                improve_conformation!(idxDprednonpred,
                                    Dij, D, P, ij_to_D, torsions,
                                    X, fixed_torsions, adj, work, num_impr, tol_lde,
                                    max_torsion
                                    )
            end

            mde, lde = MDE_LDE(Dij, D, X)

            if (lde <= tol_lde) || (mde <= tol_mde)
                stop = true
            end

            rmsdpass = true
            if !stop
                # different conformation?
                centralize!(X, Xtmp1, rmsd_idx)
                for c in 1:length(sols)
                    centralize!(sols[c], Xtmp2, rmsd_idx)
                    rmsd = rmsd_protein(Xtmp1, Xtmp2, rmsd_idx)
                    if rmsd <= tol_similar
                        rmsdpass = false
                        break
                    end
                end
                if rmsdpass
                    consec_similar = 0
                else
                    consec_similar += 1
                    if consec_similar >= max_similar
                        # maximum number of consecutive trials with similar conformations
                        stop = true
                    end
                end
            end

            spg_applied = false
            if !stop && rmsdpass
                # different infeasible conformation, try to improve it through SPG
                @views for k in 1:nd
                    i,j = Dij[k,1:2]
                    work.dists[k] = euclidean(X[1:3,i], X[1:3,j])
                end

                # apply SPG
                if verbose > 2
                    println("\n\n>>> Trying to improve conformation through SPG")
                end
                strs, spg_it, st = spg(
                    1:3*nv, idxDprednonpred, X, Dij, D, ij_to_D,
                    tol_mde, tol_lde, tol_stress,
                    max(spg_maxit, 20*nv), spg_lacktol, spg_eta, spg_lmin, spg_lmax, spg_lsm,
                    work, verbose
                )

                mde, lde = MDE_LDE(Dij, D, X)

                if (lde <= tol_lde) || (mde <= tol_mde)
                    stop = true
                end

                spg_applied = true
            end

            if rmsdpass
                push!(sols, deepcopy(X))
                push!(soltypes, -1)
                push!(ldes, lde)
                push!(mdes, mde)
                push!(stress, spg_applied ? strs : Inf)
            end

            if stop
                # X is feasible!
                soltypes[end] = spg_applied ? 2 : 1

                if verbose > 0
                    if spg_applied
                        println("\rA solution was found after applying SPG")
                    else
                        println("\rA solution was found in the initialization phase")
                    end
                end
            end
        end
        end   # end of @elapsed begin

        time_init += time_initk
        time_total = time_pre + time_init

        if stop || (time_total >= max_time)
            break
        end
    end

    if verbose > 0
        println("\r$(length(sols)) initial conformations computed in $(total_count) trials.")

        println("\nSUMMARY\n$(repeat('-',48))")
        @printf("Total time: %.6lf s\n", time_total)
        println("Number of solutions found: $(count(soltypes .>= 0))")
        b = argmin(mdes)
        @printf("Minimum MDE among all conformations: %9.3e\n", mdes[b])
        @printf("LDE of the corresponding conformation: %9.3e\n", ldes[b])
        @printf("Status of the corresponding conformation: %d\n", soltypes[b])
    end

    return sols, soltypes, ldes, mdes, time_init, time_total
end

end
