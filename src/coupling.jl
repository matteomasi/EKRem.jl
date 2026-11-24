# Coupling of transport and geochemical reactions by the sequential non-iterative (SNIA) approach:
# over each coupling step, the transport equations are integrated first, and the reactions
# are then computed in every node with the transported concentrations.

"""
The struct Coupling contains the information needed for the coupling of transport 
and geochemical reactions, i.e., a chemistry object, the columns of the transported
species in PhreeqcRM, the list of nodes where reactions are computed, and the 
PhreeqcRM selected output.
"""
struct Coupling
    chemistry::Chemistry
    columns::Vector{Int}
    nodes::Vector{Int}
    selected::Vector{Dict{String,Vector{Float64}}}
end

"""
    transported_species(chemistry, species) -> (Vector{Species}, columns)

Transported species and their columns in the PhreeqcRM species concentrations. `species` is
`nothing` (all aqueous species except H2O), a vector of names (charge numbers and diffusion
coefficients from the database), or a vector of [`Species`](@ref) (diffusion coefficients
given by the user).
"""
function transported_species(chemistry::Chemistry, species)
    all_species = aqueous_species(chemistry)
    selected = if species === nothing
        filter(s -> s.name != "H2O", all_species)
    else
        map(species) do s
            name = s isa Species ? s.name : String(s)
            i = findfirst(a -> a.name == name, all_species)
            i === nothing && throw(ArgumentError("\"$name\" is not an aqueous species of the PHREEQC system"))
            s isa Species || return all_species[i]
            s.z == all_species[i].z || throw(
                ArgumentError("The charge number of $name ($(s.z)) differs from the database ($(all_species[i].z))"),
            )
            return s
        end
    end
    columns = [findfirst(a -> a.name == s.name, all_species) for s in selected]
    return collect(Species, selected), columns
end

# Definition numbers of the initial conditions at every node: `initial` in the interior, and
# the definitions of the boundaries (reservoir or compartment) at the boundary nodes
function node_conditions(initial::NamedTuple, nnodes::Integer, left::NamedTuple, right::NamedTuple)
    for conditions in (initial, left, right), key in keys(conditions)
        key in INITIAL_KEYS ||
            throw(ArgumentError("Unknown initial condition `$key`; the keys are $(join(INITIAL_KEYS, ", "))"))
    end
    conditions = map(INITIAL_KEYS) do key
        value = get(initial, key, -1)
        value isa Integer ||
            length(value) == nnodes ||
            throw(ArgumentError("`$key` must be a number or a vector of $nnodes numbers (one per node)"))
        v = collect(Int, value isa Integer ? fill(value, nnodes) : value)
        v[1], v[end] = get(left, key, -1), get(right, key, -1)
        return v
    end
    return NamedTuple{INITIAL_KEYS}(conditions)
end

"""
    react_nodes!(coupling, state, Δt)

Compute the reactions over `Δt` in every node, with the concentrations of the transported
species in `state` (unknowns × nodes), and replace them by the concentrations after the
reactions. Negative concentrations are set to zero before they are passed to PhreeqcRM.
Boundary nodes with a reservoir keep its composition.
"""
function react_nodes!(coupling::Coupling, state::AbstractMatrix, Δt::Real)
    (; chemistry, columns, nodes) = coupling
    n = length(columns)
    c = species_concentrations(chemistry)
    c[:, columns] .= max.(transpose(state[1:n, :]), 0.0)
    set_species_concentrations!(chemistry, c)
    react!(chemistry, Δt)
    c = species_concentrations(chemistry)
    state[1:n, nodes] .= transpose(c[nodes, columns])
    return state
end

# Selected output of all cells; empty without a SELECTED_OUTPUT block in the input
function _selected_output(chemistry::Chemistry)
    LibPhreeqcRM.RM_GetSelectedOutputCount(chemistry.id) == 0 && return Dict{String,Vector{Float64}}()
    return selected_output(chemistry)
end
