using MDGP
using LinearAlgebra
using Test

@testset "Predecessor frame" begin
    X = [1.0 0.0 0.0;
         0.0 0.0 1.0;
         0.0 0.0 0.0]
    U = zeros(3, 3)

    MDGP.computeU!(U, X, 1, 2, 3)

    @test U' * U ≈ Matrix{Float64}(I, 3, 3)

    local_coordinates = [-0.25, 0.5, 0.75]
    reflected_coordinates = copy(local_coordinates)
    reflected_coordinates[3] *= -1
    candidate = X[:,1] + U * local_coordinates
    reflected_candidate = X[:,1] + U * reflected_coordinates
    candidate_distances = [norm(candidate - X[:,i]) for i in 1:3]
    reflected_distances = [norm(reflected_candidate - X[:,i]) for i in 1:3]

    @test candidate_distances ≈ reflected_distances

    predecessor_normal = cross(X[:,1] - X[:,2], X[:,3] - X[:,2])
    candidate_side = dot(predecessor_normal, candidate - X[:,2])
    reflected_side = dot(predecessor_normal, reflected_candidate - X[:,2])

    @test candidate_side ≈ -reflected_side

    collinear_X = [1.0 0.0 2.0;
                   0.0 0.0 0.0;
                   0.0 0.0 0.0]
    @test_throws ErrorException MDGP.computeU!(U, collinear_X, 1, 2, 3)

    nearly_collinear_X = [1.0 0.0 2.0;
                          0.0 0.0 1.0e-12;
                          0.0 0.0 0.0]
    @test_throws ErrorException MDGP.computeU!(U, nearly_collinear_X, 1, 2, 3)
end

@testset "Multistart (1TOS)" begin
    _, Dij, D, P, _, atoms, torsions = mdgp_read("D_1TOS.dat", "P_1TOS.dat")
    output = mdgp_multistart(Dij, D, P, atoms, torsions, verbose=0, seed=0)
    @test any(output[2] .> 0)
end
