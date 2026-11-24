# Boundary conditions at the ends of the domain (electrodes).

"Boundary condition for the concentrations and the electric potential."
abstract type BoundaryCondition end

"""
    Reservoir(concentrations; potential)
    Reservoir(concentrations; current_density)

Boundary in contact with a well-mixed reservoir of fixed composition. The concentrations
are Dirichlet conditions; `concentrations` maps species names to concentrations
(mol m⁻³), and species not listed have zero concentration. The composition must be
electroneutral.

The electric condition is one of the following:
- `potential`: fixed electric potential (V), a Dirichlet condition (voltage control);
- `current_density`: fixed electric current density (A m⁻²) applied at the boundary, per
  unit cross-sectional area of the medium (current control). It is positive at the anode,
  where the current flows into the domain, and negative at the cathode. The potential of
  the boundary is then computed. Current control requires the
  [`ChargeConservation`](@ref) formulation and a fixed potential at the other boundary.

# Examples
```julia
Reservoir(Dict("Na+" => 10.0, "Cl-" => 10.0); potential=3.0)
Reservoir(Dict("Na+" => 10.0, "Cl-" => 10.0); current_density=5.0)
```
"""
struct Reservoir <: BoundaryCondition
    concentrations::Dict{String,Float64}
    potential::Union{Float64,Nothing}
    current_density::Union{Float64,Nothing}
end

function Reservoir(
    concentrations::AbstractDict; potential::Union{Real,Nothing}=nothing, current_density::Union{Real,Nothing}=nothing
)
    xor(potential === nothing, current_density === nothing) ||
        throw(ArgumentError("A Reservoir requires either a potential or a current density, but not both"))
    return Reservoir(
        Dict(String(k) => Float64(v) for (k, v) in concentrations),
        potential === nothing ? nothing : Float64(potential),
        current_density === nothing ? nothing : Float64(current_density),
    )
end

"Whether the potential is fixed at the boundary (voltage control)."
fixed_potential(b::Reservoir) = b.potential !== nothing
