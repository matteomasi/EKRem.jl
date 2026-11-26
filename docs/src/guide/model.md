# Model set-up

A one-dimensional model is defined by [`Model1D`](@ref) with keyword arguments: the
geometry and the grid, the transported species, the initial composition, the boundary
conditions, the formulation of the potential and the properties of the medium. All
quantities are in SI units; concentrations are in mol per m³ of pore water.

## Grid

`grid` is either a vector of node coordinates from 0 to `length`, or a named tuple for a
grid refined toward both ends ([`grid_refined`](@ref)):

```julia
grid = collect(range(0, 0.1; length=201))          # uniform, 0.5 mm
grid = (h_min=1.0e-5, h_max=1.0e-3, growth=1.1)    # refined toward the boundaries
```

Sharp gradients form near the boundaries, at the reservoirs and the compartments, so
refined grids are usually more efficient than uniform ones.

## Species and composition

The transported species are given as a vector of [`Species`](@ref), with name, charge
number and diffusion coefficient in free water (m² s⁻¹). The initial composition is uniform
and given as a dictionary of concentrations; species that are not listed are absent. All
compositions must be electroneutral.

```julia
species = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9)]
initial = Dict("Na+" => 10.0, "Cl-" => 10.0)       # mol/m³
```

In a model with geochemical reactions, the species and the compositions are taken from the
PHREEQC database and input (see [Geochemical reactions](@ref reactions-guide)).

## Porous medium

`porosity`, `tortuosity`, `relative_permittivity` and `electroosmotic_permeability` define
the medium ([`PorousMedium`](@ref)). The effective diffusion coefficient of each species is
``D_i^* = n \tau D_i``. With the default values (``n = \tau = 1``, no electroosmosis), the
equations describe a free solution. `temperature` (K) enters the ratio ``F/(R T)`` of the
electromigration term.

With a non-zero `electroosmotic_permeability` ``k_{eo}`` (m² V⁻¹ s⁻¹), the pore water flows
with the uniform flux ``q = k_{eo} [\phi(0) - \phi(L)] / L``, from the anode to the cathode
for ``k_{eo} > 0`` (see [Electroosmosis](@ref electroosmosis)). The flux is returned by
[`electroosmotic_flow`](@ref).

## Boundary conditions

Each end of the column is a [`Reservoir`](@ref) of fixed composition or an electrode
[`Compartment`](@ref). The electric condition is a fixed potential (`potential`, V) or an
applied current density (`current_density`, A m⁻²):

| Control | Left boundary | Right boundary | Formulation of the potential |
|:--------|:--------------|:---------------|:-----------------------------|
| Voltage | `potential` | `potential` | any |
| Current | `current_density` | `potential` | [`ChargeConservation`](@ref) |
| Current | `potential` | `current_density` | [`ChargeConservation`](@ref) |

The current density is positive where the current enters the column (anode) and negative
where it leaves it (cathode).

Electrode compartments require the cross-sectional `area` of the column, which converts the
compartment volume and the flushing rate to quantities per unit area. They can be used under
current or voltage control. The H⁺ and OH⁻ produced by electrolysis must be among the
transported species (both of them under voltage control, where the anode follows from the
direction of the computed current).

```julia
model = Model1D(;
    length=0.1,
    grid=(h_min=1.0e-5, h_max=1.0e-3),
    species=[Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9),
             Species("H+", 1, 9.31e-9), Species("OH-", -1, 5.27e-9)],
    initial=Dict("Na+" => 10.0, "Cl-" => 10.0),
    left=Compartment(Dict("Na+" => 10.0, "Cl-" => 10.0); volume=1.0e-4, current_density=10.0),
    right=Compartment(Dict("Na+" => 10.0, "Cl-" => 10.0); volume=1.0e-4, potential=0.0),
    potential=ChargeConservation(),
    porosity=0.4,
    tortuosity=0.5,
    area=1.0e-3,
)
```

## Formulation of the potential

[`Poisson`](@ref) (default) resolves the charge separation in boundary layers whose
thickness is of the order of the Debye length. When the Debye length is much smaller than
the smallest grid spacing ([`debye_ratio`](@ref) much smaller than 1), the layers are not
resolved, and [`ChargeConservation`](@ref) gives the same electroneutral solution; it is
also required under current control. [`FixedField`](@ref) decouples the
potential from the concentrations.
