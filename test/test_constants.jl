using Test
using EKRem

@testset "Constants and electrolyte scales" begin
    @test FARADAY ≈ 96485.33212 rtol = 1.0e-10
    @test GAS_CONSTANT ≈ 8.314462618 rtol = 1.0e-9

    # Debye length of a 10 mmol/L 1:1 electrolyte in water at 25 °C: about 3.04 nm, 
    # using the approximation λ_D ≈ 0.304 nm / √(I / (mol/L)).
    @test debye_length([10.0, 10.0], [1, -1]) ≈ 3.04e-9 rtol = 0.01
    # λ_D scales with √εr
    @test debye_length([10.0, 10.0], [1, -1]; εr=1.0) ≈
        debye_length([10.0, 10.0], [1, -1]) / sqrt(WATER_RELATIVE_PERMITTIVITY)

    # Charge relaxation time ε/σ of 0.1 mol/L NaCl: σ ≈ 1.07 S/m (Nernst–Einstein estimate),
    # τ ≈ 0.65 ns
    τ = charge_relaxation_time([100.0, 100.0], [1, -1], [1.33e-9, 2.03e-9])
    @test 5.0e-10 < τ < 8.0e-10
end
