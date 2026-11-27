using Test
using SimEK

@testset "Physical constants" begin
    @test FARADAY ≈ 96485.33212 rtol = 1.0e-10
    @test GAS_CONSTANT ≈ 8.314462618 rtol = 1.0e-9
end
