# Properties of electrolyte solutions: length and time scales of the electric charge, and
# conductivity.

"""
    debye_length(c, z; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)

Debye length (m) of an aqueous solution with concentrations `c` (mol m⁻³) and charge
numbers `z`:

    λ_D = √(ε₀ εr R T / (F² Σ zᵢ² cᵢ))

`εr` is the relative permittivity of the pore water. In the bulk form of the Poisson
equation, `εr` corresponds to ε_eff / (n ε₀).
"""
function debye_length(c, z; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)
    s = sum(zi^2 * ci for (zi, ci) in zip(z, c))
    return sqrt(VACUUM_PERMITTIVITY * εr * GAS_CONSTANT * T / (FARADAY^2 * s))
end

"""
    conductivity(c, z, D; T=298.15)

Electric conductivity (S m⁻¹) of an aqueous solution with concentrations `c` (mol m⁻³),
charge numbers `z` and diffusion coefficients `D` (m² s⁻¹), from the Nernst–Einstein
relation:

    σ = F² / (R T) Σ zᵢ² Dᵢ cᵢ

With the effective diffusion coefficients D* = n τ D of a porous medium, the result is the
bulk conductivity of the water-saturated medium.
"""
function conductivity(c, z, D; T=298.15)
    return FARADAY^2 / (GAS_CONSTANT * T) * sum(zi^2 * Di * ci for (zi, Di, ci) in zip(z, D, c))
end

"""
    charge_relaxation_time(c, z, D; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)

Charge (dielectric) relaxation time τ = ε / σ (s) of an aqueous solution, with the
conductivity σ of [`conductivity`](@ref). Concentrations `c` in mol m⁻³, diffusion
coefficients `D` in m² s⁻¹.
"""
function charge_relaxation_time(c, z, D; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)
    return VACUUM_PERMITTIVITY * εr / conductivity(c, z, D; T)
end
