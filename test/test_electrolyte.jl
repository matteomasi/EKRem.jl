using Test
using SimEK

@testset "Properties of electrolyte solutions" begin
    # Debye length of a 10 mmol/L 1:1 electrolyte in water at 25 °C: about 3.04 nm, 
    # using the approximation λ_D ≈ 0.304 nm / √(I / (mol/L)).
    @test debye_length([10.0, 10.0], [1, -1]) ≈ 3.04e-9 rtol = 0.01
    # λ_D scales with √εr
    @test debye_length([10.0, 10.0], [1, -1]; εr=1.0) ≈
        debye_length([10.0, 10.0], [1, -1]) / sqrt(WATER_RELATIVE_PERMITTIVITY)

    # Nernst–Einstein: molar conductivity λ = F² D / (R T) = 50.1 S cm² mol⁻¹ for D = 1.334e-9 m²/s,
    # i.e. σ = 5.01e-3 S/m at 1 mol/m³; σ is proportional to z² and to the concentration
    @test conductivity([1.0], [1], [1.334e-9]) ≈ 5.01e-3 rtol = 1.0e-3
    @test conductivity([3.0], [2], [1.0e-9]) ≈ 12 * conductivity([1.0], [1], [1.0e-9])

    # Charge relaxation time ε/σ of 0.1 mol/L NaCl: σ ≈ 1.26 S/m (Nernst–Einstein), τ ≈ 0.55 ns
    τ = charge_relaxation_time([100.0, 100.0], [1, -1], [1.33e-9, 2.03e-9])
    @test τ ≈
        VACUUM_PERMITTIVITY * WATER_RELATIVE_PERMITTIVITY / conductivity([100.0, 100.0], [1, -1], [1.33e-9, 2.03e-9])
    @test 5.0e-10 < τ < 6.0e-10
end
