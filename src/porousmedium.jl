# Properties of the water-saturated porous medium.

@doc raw"""
    PorousMedium(; porosity=1.0, tortuosity=1.0, relative_permittivity=WATER_RELATIVE_PERMITTIVITY,
                 electroosmotic_permeability=0.0)

Properties of a water-saturated porous medium:

- `porosity` ``n`` (–);
- `tortuosity` factor ``\tau`` (–), with ``0 < \tau \le 1``, so that the effective
  diffusion coefficient of a species is ``D^* = n \tau D``;
- `relative_permittivity` ``\varepsilon_r`` of the pore water (–). The effective
  permittivity of the medium in the Poisson equation is
  ``\varepsilon_\mathrm{eff} = n \varepsilon_0 \varepsilon_r``;
- `electroosmotic_permeability` ``k_{eo}`` (m² V⁻¹ s⁻¹), per unit cross-sectional area of
  the medium, so that the electroosmotic water flux (Darcy flux) is
  ``q = -k_{eo} \nabla\phi``. It is positive when the flow is directed from the anode to the
  cathode, as for solids with a negative zeta potential ``\zeta``. By the
  Helmholtz–Smoluchowski equation with the correction of Casagrande,
  ``k_{eo} = -\varepsilon \zeta n \tau / \mu`` (Masi et al., 2017, Eq. 6).

With the default values (``n = \tau = 1``, ``k_{eo} = 0``) the equations describe a free
solution.
"""
struct PorousMedium
    porosity::Float64
    tortuosity::Float64
    relative_permittivity::Float64
    electroosmotic_permeability::Float64
    function PorousMedium(;
        porosity::Real=1.0,
        tortuosity::Real=1.0,
        relative_permittivity::Real=WATER_RELATIVE_PERMITTIVITY,
        electroosmotic_permeability::Real=0.0,
    )
        0 < porosity <= 1 || throw(ArgumentError("The porosity must lie in (0, 1], got $porosity"))
        0 < tortuosity <= 1 || throw(ArgumentError("The tortuosity factor must lie in (0, 1], got $tortuosity"))
        relative_permittivity > 0 ||
            throw(ArgumentError("The relative permittivity must be positive, got $relative_permittivity"))
        isfinite(electroosmotic_permeability) ||
            throw(ArgumentError("The electroosmotic permeability must be finite, got $electroosmotic_permeability"))
        return new(porosity, tortuosity, relative_permittivity, electroosmotic_permeability)
    end
end

@doc raw"Effective permittivity ``\varepsilon_\mathrm{eff} = n \varepsilon_0 \varepsilon_r`` (F m⁻¹) of the medium."
effective_permittivity(m::PorousMedium) = m.porosity * VACUUM_PERMITTIVITY * m.relative_permittivity

@doc raw"Effective diffusion coefficient ``D^* = n \tau D`` (m² s⁻¹) of a species in the medium."
effective_diffusion(m::PorousMedium, s::Species) = m.porosity * m.tortuosity * s.D
