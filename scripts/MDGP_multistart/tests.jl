###########################
# Simple script to reproduce the tests with Algorithm 3 as reported in
# Secchin, da Rocha, da Rosa, Liberti, Lavor. A fast heuristic for the molecular distance geometry problem. 2025
###########################

using MDGP
using Base.Threads
using Printf
using LinearAlgebra
using DataFrames

BLAS.set_num_threads(1)

function run_multistart()
    # dataset folder
    dataset = "dataset"

    if !isdir("output")
        mkdir("output")
    end
    if !isdir("data")
        mkdir("data")
    end
    if !isdir(dataset)
        mkdir(dataset)
    end

    # generate instances from PDB (this is done only once)
    for p in readlines("selected_probs")
        Dfile = "$(dataset)/$(p)/I_$(p)_model1_chainA_ddgpHCorder1.dat"
        Pfile = "$(dataset)/$(p)/T_$(p)_model1_chainA_ddgpHCorder1.dat"
        Xfile = "$(dataset)/$(p)/X_$(p)_model1_chainA_ddgpHCorder1.dat"

        if !isfile(Dfile) || !isfile(Pfile) || !isfile(Xfile)
            println("\nParsing $(p)\n$(repeat('=',50))\n")
            try
                run(`python3 pdb_parser.py --data_dir data --output_dir output --pdb_id $p --model 1 --chain A --ddgp_hc_order 1 --cut 5.0 --interval_width 2.0 --local_interval_width 1.0 --angular_width 50 1`)
            catch
                println("Error while parsing $(p)")
            end
        end
    end

    # delete temporary files
    rm("data", force=true, recursive=true)
    rm("output", force=true, recursive=true)

    # read generated instances
    probs = readdir(dataset)

    @assert !isempty(probs) "No problems to solve! Verify that the PDB parser is working properly"

    results = DataFrame(
                ID=String[],
                n=Int64[],
                E=Int64[],
                Nc=Int64[],
                status=Int64[],
                LDE=Float64[],
                MDE=Float64[],
                time=Float64[]
                )

    applock = SpinLock()

    Threads.@threads for p in probs
        println("\nSolving $(p)\n$(repeat('=',50))\n")
        try
            Dfile = "$(dataset)/$(p)/I_$(p)_model1_chainA_ddgpHCorder1.dat"
            Pfile = "$(dataset)/$(p)/T_$(p)_model1_chainA_ddgpHCorder1.dat"
            Xfile = "$(dataset)/$(p)/X_$(p)_model1_chainA_ddgpHCorder1.dat"
            _, Dij, D, P, _, atoms, torsions = MDGP_read(Dfile, Pfile, Xfile=Xfile)

            # remove sign in P for "quasi-planar" atoms
            planarity!(Dij, P, torsions)

            # solve p
            _, status, ldes, mdes, _, time =
                MDGP_multistart(Dij, D, P, atoms, torsions, verbose=0)

            # best conformation
            b = argmin(mdes)

            # save partial results
            lock(applock)
            push!(results, (p, size(P,1), size(D,1), length(ldes), status[b], ldes[b], mdes[b], time))
            unlock(applock)
        catch
            println("Error solving $(p)")
        end
    end

    # sort results
    sort!(results, [2,3,1])

    # export results
    write("RESULTS", @sprintf("%s", results))
end

# if the reference torsion angle plus/minus 'tol_planar'
# contains 0 or 180 degrees, zeroing the correspondent sign in P
function planarity!(Dij, P, torsions; tol_planar = 20.0)
    nv = size(P,1)
    nd = size(Dij,1)

    @inbounds @views for i in 4:nv
        if P[i,4] != 0
            torsion, delta = torsions[i,1:2]

            if delta == 0.0
                delta = tol_planar
            end

            if ((torsion - delta < 0.0) && (torsion + delta > 0.0)) ||
                ((torsion - delta < 180.0) && (torsion + delta > 180.0))
                P[i,4] = 0
            end
        end
    end
end
