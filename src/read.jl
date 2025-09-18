# compute torsions angles from X
function omega(l, P, X)
    l1, l2, l3 = P[l,1:3]

    a = X[1:3,l2] - X[1:3,l3]
    b = X[1:3,l1] - X[1:3,l2]
    c = X[1:3,l ] - X[1:3,l1]

    return 180 * atan(
        norm(b) * dot(a, cross(b, c)),
        dot(cross(a, b), cross(b, c))
        ) / pi
end

"""
    X, Dij, D, P, res, atoms, torsions = MDGP_read(dir; fileid=1, model=1, chain='A')
"""
function MDGP_read(
    dir::String;
    fileid = 1,
    model = 1,
    chain = 'A'
)
    if !isdir(dir)
        return [], [], [], [], [], [], []
    end

    files = readdir(dir)

    Dij = []
    D = []
    X = []
    P = []
    res = []
    atoms = []
    resnames = []
    angles = []

    iI = 0
    iX = 0
    iT = 0
    if fileid < 11
        for i in 1:length(files)
            if contains(files[i], "model$(model)_chain$(chain)") && contains(files[i], "ddgpHCorder$(fileid).dat")
                if startswith(files[i], "I_")
                    iI = i
                elseif startswith(files[i], "X_")
                    iX = i
                elseif startswith(files[i], "T_")
                    iT = i
                end
            end
        end
    else
        for i in 1:length(files)
            if contains(files[i], "model$(model)_chain$(chain)") && contains(files[i], "dmdgpHCorder.dat")
                if startswith(files[i], "I_")
                    iI = i
                elseif startswith(files[i], "X_")
                    iX = i
                elseif startswith(files[i], "T_")
                    iT = i
                end
            end
        end
    end

    if iI > 0
        dadosI = readdlm(joinpath(dir,files[iI]))

        # read columns of "I file"
        D = Float64.([dadosI[:,5] dadosI[:,6]])

        # Obs: the second indice is read first
        Dij = Int64.([dadosI[:,2] dadosI[:,1]])
        res = Int64.([dadosI[:,4] dadosI[:,3]])
        atoms = [dadosI[:,8] dadosI[:,7]]
        resnames = [dadosI[:,10] dadosI[:,9]]
    else
        @error "File 'I' not found."
        return [], [], [], [], [], [], []
    end

    n = maximum(Dij)

    if iX > 0
        dadosX = readdlm(joinpath(dir,files[iX]))

        X = Float64.(dadosX[:,1:3])
    else
        @warn "File 'X' not found."
    end

    P = Matrix{Int64}(undef, n, 4)

    if iT > 0
        dadosT = readdlm(joinpath(dir,files[iT]))
        P[:,1:4] = Int64.(dadosT[:,2:5])
        # recompute signs if X is available
        if iX > 0
            @inbounds for i in 4:n
                if P[i,4] != 0
                    P[i,4] = Int64(sign(omega(i,P,X')))
                end
            end
        end
        torsion_angles = abs.(Float64.(dadosT[:,6:7]))
    else
        @warn "File 'T' not found."

        # P was not provided. Consider the DMDGP order
        P = Matrix{Int64}(undef, n, 4)
        @views for i in 1:n
            P[i,1:4] .= [max(i-1,0); max(i-2,0); max(i-3,0); 0]
        end
        torsion_angles = zeros(n, 2)
    end

    idxDpred = Int64[]      # between predecessors
    idxDnonpred = Int64[]   # extra
    idxDvdwc = Int64[]      # Van der Waals between consecutive residues
    idxDvdwnc = Int64[]     # other Van der Waals
    @inbounds for k in 1:size(D,1)
        i,j = Dij[k,1:2]
        if (j in P[i,1:3]) || (i in P[j,1:3])
            push!(idxDpred, k)
        elseif D[k,2] < 900.0
            push!(idxDnonpred, k)
        elseif abs(res[k,1] - res[k,2]) <= 1
            push!(idxDvdwc, k)
        else
            push!(idxDvdwnc, k)
        end
    end

    # sort distances
    TMP = sortslices([Dij D res atoms resnames][idxDpred,:], dims=1, by=row -> (row[2],row[1]))
    TMP = [TMP; sortslices([Dij D res atoms resnames][idxDnonpred,:], dims=1, by=row -> (row[2],row[1]))]
    TMP = [TMP; sortslices([Dij D res atoms resnames][idxDvdwc,:], dims=1, by=row -> (row[2],row[1]))]
    TMP = [TMP; sortslices([Dij D res atoms resnames][idxDvdwnc,:], dims=1, by=row -> (row[2],row[1]))]

    Dij = Int64.(TMP[:,1:2])
    D   = Float64.(TMP[:,3:4])
    res = Int64.(TMP[:,5:6])
    atoms = TMP[:,7:8]
    resnames = TMP[:,9:10]

    # make explicit infinite upper bounds
    D[ D .>= 900.0 ] .= Inf

    # map between atoms indices to their types and residues
    atoms_map = Vector{String}(undef, n)
    res_map = zeros(Int64, n)
    atoms_map .= "?"
    @inbounds @views for k in 1:size(Dij,1)
        atoms_map[Dij[k,1]] = strip(atoms[k,1])
        atoms_map[Dij[k,2]] = strip(atoms[k,2])
        res_map[Dij[k,1]] = res[k,1]
        res_map[Dij[k,2]] = res[k,2]
    end

    return Matrix(X'), Dij, D, P, res_map, atoms_map, torsion_angles
end
