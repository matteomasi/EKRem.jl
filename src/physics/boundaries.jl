# Boundary conditions at the ends of the domain (electrodes).

"Boundary condition for the concentrations and the electric potential."
abstract type BoundaryCondition end

"""
    Reservoir(concentrations; potential)
    Reservoir(concentrations; current_density)
    Reservoir(; solution, potential)
    Reservoir(; solution, current_density)

Boundary in contact with a well-mixed reservoir of fixed composition. The concentrations
are Dirichlet conditions; `concentrations` maps species names to concentrations
(mol m⁻³), and species not listed have zero concentration. The composition must be
electroneutral. In a model with geochemical reactions, the composition is instead given by
the number of a `solution` defined in the PHREEQC input.

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
Reservoir(; solution=2, potential=3.0)
```
"""
struct Reservoir <: BoundaryCondition
    concentrations::Dict{String,Float64}
    solution::Union{Int,Nothing}
    potential::Union{Float64,Nothing}
    current_density::Union{Float64,Nothing}
end

function Reservoir(
    concentrations::AbstractDict; potential::Union{Real,Nothing}=nothing, current_density::Union{Real,Nothing}=nothing
)
    return _reservoir(Dict(String(k) => Float64(v) for (k, v) in concentrations), nothing, potential, current_density)
end

function Reservoir(;
    solution::Integer, potential::Union{Real,Nothing}=nothing, current_density::Union{Real,Nothing}=nothing
)
    return _reservoir(Dict{String,Float64}(), Int(solution), potential, current_density)
end

function _reservoir(concentrations, solution, potential, current_density)
    return Reservoir(concentrations, solution, _electric_condition("Reservoir", potential, current_density)...)
end

function _electric_condition(type, potential, current_density)
    xor(potential === nothing, current_density === nothing) ||
        throw(ArgumentError("A $type requires either a potential or a current density, but not both"))
    return potential === nothing ? nothing : Float64(potential),
    current_density === nothing ? nothing : Float64(current_density)
end

@doc raw"""
    Compartment(concentrations; volume, flow=0.0, current_density)
    Compartment(concentrations; volume, flow=0.0, potential)
    Compartment(; initial, volume, flow=0.0, current_density)
    Compartment(; initial, volume, flow=0.0, potential)

Electrode compartment at the boundary: a well-mixed solution of `volume` ``V`` (m³) in
contact with the end of the column, flushed with deionised water at the rate `flow` ``Q``
(m³ s⁻¹; the same volume flows out). The concentrations in the compartment change by the
fluxes from and into the column, by the outflow and by water electrolysis at the electrode:

```math
V \frac{\mathrm{d} c_i}{\mathrm{d} t} = \mp A J_i - Q c_i + R_i,
```

with ``A`` the cross-sectional area of the column (the `area` of the model) and ``J_i`` the
flux at the boundary in the ``+x`` direction (``-`` at the left, ``+`` at the right
boundary). Electrolysis produces H⁺ at the anode
(``\mathrm{2\,H_2O \to O_2 + 4\,H^+ + 4\,e^-}``) and OH⁻ at the cathode
(``\mathrm{2\,H_2O + 2\,e^- \to H_2 + 2\,OH^-}``), at the rate ``R_i = I / F``, where ``I``
is the current; the gases are assumed to escape.

The initial composition is given by `concentrations` (mol m⁻³) or, in a model with
geochemical reactions, by the numbers of PHREEQC definitions, e.g.
`initial=(solution=2, equilibrium_phases=2)`; with geochemical reactions, the solution of
the compartment reacts at every coupling step, like the pore water.

The electric condition is given as for a [`Reservoir`](@ref). Under current control, the
electrolysis rate follows from the applied current, which is the same at both electrodes.
Under voltage control (potential fixed at both boundaries), it follows from the computed
current between the compartment and the column, so that the compartment remains
electroneutral; the anode is then the electrode where the current enters the column, and
both H⁺ and OH⁻ must be species of the model.
"""
struct Compartment <: BoundaryCondition
    concentrations::Dict{String,Float64}
    initial::Union{NamedTuple,Nothing}
    volume::Float64
    flow::Float64
    potential::Union{Float64,Nothing}
    current_density::Union{Float64,Nothing}
end

function Compartment(
    concentrations::AbstractDict;
    volume::Real,
    flow::Real=0.0,
    potential::Union{Real,Nothing}=nothing,
    current_density::Union{Real,Nothing}=nothing,
)
    c = Dict(String(k) => Float64(v) for (k, v) in concentrations)
    return _compartment(c, nothing, volume, flow, potential, current_density)
end

function Compartment(;
    initial::NamedTuple,
    volume::Real,
    flow::Real=0.0,
    potential::Union{Real,Nothing}=nothing,
    current_density::Union{Real,Nothing}=nothing,
)
    return _compartment(Dict{String,Float64}(), initial, volume, flow, potential, current_density)
end

function _compartment(concentrations, initial, volume, flow, potential, current_density)
    volume > 0 || throw(ArgumentError("The volume of a compartment must be positive"))
    flow >= 0 || throw(ArgumentError("The flushing rate of a compartment must be non-negative"))
    electric = _electric_condition("Compartment", potential, current_density)
    return Compartment(concentrations, initial, Float64(volume), Float64(flow), electric...)
end

"Whether the potential is fixed at the boundary (voltage control)."
fixed_potential(b::BoundaryCondition) = b.potential !== nothing

# PHREEQC definitions of the boundary node in a model with geochemical reactions
function phreeqc_definitions(b::Reservoir)
    b.solution === nothing && throw(ArgumentError("With geochemical reactions, reservoirs are PHREEQC solutions"))
    return (solution=b.solution,)
end

function phreeqc_definitions(b::Compartment)
    b.initial === nothing &&
        throw(ArgumentError("With geochemical reactions, compartments are given by PHREEQC definitions (`initial`)"))
    return b.initial
end
