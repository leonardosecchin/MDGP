using MDGP
using Test

@testset "1TOS" begin
    _, Dij, D, P, res, atoms, torsions = MDGP_read(".", fileid=1)
    output = MDGP_multistart(Dij, D, P, atoms, res, torsions, verbose=0)
    @test any(output[2] .> 0)
end
