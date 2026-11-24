using Test
using EKRem

@testset "Grids refined toward both ends" begin
    L = 0.01
    x = grid_refined(L; h_min=1.0e-5, h_max=4.0e-5, growth=1.1)
    h = diff(x)
    @test x[1] == 0 && x[end] == L
    @test all(h .> 0)
    # Symmetric with respect to the centre, which is a node
    @test x ≈ L .- reverse(x)
    @test any(isapprox.(x, L / 2))
    # Smallest cells at the ends; regular cells between h_min and h_max, growing by at most
    # the growth factor (the adjusted central cells excepted)
    @test h[1] ≈ 1.0e-5 && h[end] ≈ 1.0e-5
    regular = h[2:(end ÷ 2 - 1)]
    @test all(1.0e-5 .<= regular .<= 4.0e-5 * (1 + 1.0e-12))
    @test all(regular[2:end] ./ regular[1:(end - 1)] .<= 1.1 + 1.0e-12)
    # Uniform grid when h_max = h_min
    u = grid_refined(1.0; h_min=0.1)
    @test u ≈ 0:0.1:1
    # Invalid arguments
    @test_throws ArgumentError grid_refined(L; h_min=0.01)
    @test_throws ArgumentError grid_refined(L; h_min=1.0e-5, h_max=1.0e-6)
end
