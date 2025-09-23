using MDGP
using Test

@testset "Multistart (1TOS)" begin
    _, Dij, D, P, _, atoms, torsions = mdgp_read("D_1TOS.dat", "P_1TOS.dat")
    output = mdgp_multistart(Dij, D, P, atoms, torsions, verbose=0, seed=0)
    @test any(output[2] .> 0)
end
