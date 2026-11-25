# Properties of electrolyte solutions: length and time scales of the electric charge, and
# conductivity.

@doc raw"""
    debye_length(c, z; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)

Debye length (m) of an aqueous solution with concentrations `c` (mol m⁻³) and charge
numbers `z`:

```math
\lambda_D = \sqrt{\frac{\varepsilon_0 \varepsilon_r R T}{F^2 \sum_i z_i^2 c_i}}
```

`εr` is the relative permittivity ``\varepsilon_r`` of the pore water. In the bulk form of
the Poisson equation, ``\varepsilon_r`` corresponds to ``\varepsilon_\mathrm{eff} / (n \varepsilon_0)``.
"""
function debye_length(c, z; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)
    s = sum(zi^2 * ci for (zi, ci) in zip(z, c))
    return sqrt(VACUUM_PERMITTIVITY * εr * GAS_CONSTANT * T / (FARADAY^2 * s))
end

@doc raw"""
    conductivity(c, z, D; T=298.15)

Electric conductivity (S m⁻¹) of an aqueous solution with concentrations `c` (mol m⁻³),
charge numbers `z` and diffusion coefficients `D` (m² s⁻¹), from the Nernst–Einstein
relation:

```math
\sigma = \frac{F^2}{R T} \sum_i z_i^2 D_i c_i
```

With the effective diffusion coefficients ``D_i^* = n \tau D_i`` of a porous medium, the
result is the bulk conductivity of the water-saturated medium.
"""
function conductivity(c, z, D; T=298.15)
    return FARADAY^2 / (GAS_CONSTANT * T) * sum(zi^2 * Di * ci for (zi, Di, ci) in zip(z, D, c))
end

@doc raw"""
    charge_relaxation_time(c, z, D; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)

Charge (dielectric) relaxation time ``\varepsilon_0 \varepsilon_r / \sigma`` (s) of an
aqueous solution, with the conductivity ``\sigma`` of [`conductivity`](@ref). Concentrations `c` in mol m⁻³, diffusion
coefficients `D` in m² s⁻¹.
"""
function charge_relaxation_time(c, z, D; εr=WATER_RELATIVE_PERMITTIVITY, T=298.15)
    return VACUUM_PERMITTIVITY * εr / conductivity(c, z, D; T)
end
