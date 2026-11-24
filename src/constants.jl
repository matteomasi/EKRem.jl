# Physical constants (SI units) and electrolyte length and time scales.

const AVOGADRO = 6.02214076e23              # mol⁻¹
const ELEMENTARY_CHARGE = 1.602176634e-19   # C
const BOLTZMANN = 1.380649e-23              # J K⁻¹

"Faraday constant F = N_A e (C mol⁻¹)."
const FARADAY = AVOGADRO * ELEMENTARY_CHARGE

"Molar gas constant R = N_A k (J mol⁻¹ K⁻¹)."
const GAS_CONSTANT = AVOGADRO * BOLTZMANN

"Vacuum permittivity ε₀ (F m⁻¹), CODATA 2018 recommended value."
const VACUUM_PERMITTIVITY = 8.8541878128e-12

"""
Relative permittivity of water at 25 °C (78.30), measured by Malmberg and Maryott (1956),
J. Res. Natl. Bur. Stand. 56(1). The temperature dependence is not implemented.
"""
const WATER_RELATIVE_PERMITTIVITY = 78.30

_ionic_sum(f, c, z) = sum(f(zi) * ci for (zi, ci) in zip(z, c))

"""
    debye_length(c, z; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)

Debye length (m) of an aqueous solution with concentrations `c` (mol m⁻³) and charge
numbers `z`:

    λ_D = √(ε₀ εr R T / (F² Σ zᵢ² cᵢ))

`εr` is the relative permittivity of the pore water. In the bulk form of the Poisson
equation, `εr` corresponds to ε_eff / (n ε₀).
"""
function debye_length(c, z; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)
    s = _ionic_sum(zi -> zi^2, c, z)
    return sqrt(VACUUM_PERMITTIVITY * εr * GAS_CONSTANT * T / (FARADAY^2 * s))
end

"""
    charge_relaxation_time(c, z, D; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)

Charge (dielectric) relaxation time τ = ε / σ (s) of an aqueous solution, with the
conductivity σ = F²/(R T) Σ zᵢ² Dᵢ cᵢ (S m⁻¹) obtained from the Nernst–Einstein relation.
Concentrations `c` in mol m⁻³, diffusion coefficients `D` in m² s⁻¹.
"""
function charge_relaxation_time(c, z, D; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)
    σ = FARADAY^2 / (GAS_CONSTANT * T) * sum(zi^2 * Di * ci for (zi, Di, ci) in zip(z, D, c))
    return VACUUM_PERMITTIVITY * εr / σ
end
