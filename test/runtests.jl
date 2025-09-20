using MDGP
using Test

@testset "1TOS" begin
    _, Dij, D, P, res, atoms, torsions = MDGP_read("D_1TOS.dat", "P_1TOS.dat")
    output = MDGP_multistart(Dij, D, P, atoms, res, torsions, verbose=0)
    @test any(output[2] .> 0)
end
