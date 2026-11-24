# Properties of the water-saturated porous medium.

"""
    PorousMedium(; porosity=1.0, tortuosity=1.0, relative_permittivity=WATER_RELATIVE_PERMITTIVITY)

Properties of a water-saturated porous medium:

- `porosity` n (–);
- `tortuosity` factor τ (–), with 0 < τ ≤ 1, so that the effective diffusion coefficient
  of a species is D* = n τ D;
- `relative_permittivity` εr of the pore water (–). The effective permittivity of the
  medium in the Poisson equation is ε_eff = n ε₀ εr.

With the default values (n = τ = 1) the equations describe a free solution.
"""
struct PorousMedium
    porosity::Float64
    tortuosity::Float64
    relative_permittivity::Float64
    function PorousMedium(;
        porosity::Real=1.0, tortuosity::Real=1.0, relative_permittivity::Real=WATER_RELATIVE_PERMITTIVITY
    )
        0 < porosity <= 1 || throw(ArgumentError("The porosity must lie in (0, 1], got $porosity"))
        0 < tortuosity <= 1 || throw(ArgumentError("The tortuosity factor must lie in (0, 1], got $tortuosity"))
        relative_permittivity > 0 ||
            throw(ArgumentError("The relative permittivity must be positive, got $relative_permittivity"))
        return new(porosity, tortuosity, relative_permittivity)
    end
end

"Effective permittivity ε_eff = n ε₀ εr (F m⁻¹) of the medium."
effective_permittivity(m::PorousMedium) = m.porosity * VACUUM_PERMITTIVITY * m.relative_permittivity

"Effective diffusion coefficient D* = n τ D (m² s⁻¹) of a species in the medium."
effective_diffusion(m::PorousMedium, s::Species) = m.porosity * m.tortuosity * s.D
