# Boundary conditions at the ends of the domain (electrodes).

"Boundary condition for the concentrations and the electric potential."
abstract type BoundaryCondition end

"""
    Reservoir(concentrations; potential)

Boundary in contact with a well-mixed reservoir of fixed composition, kept at a fixed
electric potential (V): the concentrations and the potential are Dirichlet conditions. 
`concentrations` maps species names to concentrations (mol m⁻³); species not
listed have zero concentration. The composition must be electroneutral.

# Example
```julia
Reservoir(Dict("Na+" => 10.0, "Cl-" => 10.0); potential=3.0)
```
"""
struct Reservoir <: BoundaryCondition
    concentrations::Dict{String,Float64}
    potential::Float64
end

function Reservoir(concentrations::AbstractDict; potential::Real)
    return Reservoir(Dict(String(k) => Float64(v) for (k, v) in concentrations), Float64(potential))
end
