# Geochemical reactions computed by PhreeqcRM: initial conditions, transfer of
# concentrations, and reaction steps.

const INITIAL_KEYS = (:solution, :equilibrium_phases, :exchange, :surface, :gas_phase, :solid_solutions, :kinetics)

"""
    Chemistry(ncells; database, input, initial, threads=1, output=nothing)

Geochemical reactions in `ncells` cells, computed by PhreeqcRM. Each cell contains 1 L of
pore water, and the amounts of minerals and other reactants defined in the PHREEQC input
are taken per litre of pore water (not per kilogram of water as in PHREEQC). Concentrations
are exchanged as concentrations of the aqueous species, in mol per m³ of pore water.

# Arguments
- `database`: PHREEQC database file.
- `input`: PHREEQC input with the definitions of the solutions, reactants and selected
  output, given as the path of a file or as text.
- `initial`: definition numbers of the initial conditions, as a named tuple with keys among
  `solution`, `equilibrium_phases`, `exchange`, `surface`, `gas_phase`, `solid_solutions`
  and `kinetics`, e.g. `(solution=1, equilibrium_phases=1)`. Each value is a number for
  all cells or a vector with one number per cell; omitted entities are absent.
- `threads`: number of threads (≤ 0 for all processors).
- `output`: path and prefix of the PhreeqcRM output files (`.chem.txt`, `.log.txt`). No
  files are written by default.

The cells are equilibrated when the object is created. The PhreeqcRM instance is released
with `close(chemistry)`, or when the object is garbage-collected.
"""
mutable struct Chemistry
    id::Int
    ncells::Int
    components::Vector{String}
    species::Vector{String}
    time::Float64
    function Chemistry(id::Integer, ncells::Integer)
        chemistry = new(id, ncells, String[], String[], 0.0)
        return finalizer(close, chemistry)
    end
end

function Chemistry(
    ncells::Integer;
    database::AbstractString,
    input::AbstractString,
    initial::NamedTuple,
    threads::Integer=1,
    output::Union{AbstractString,Nothing}=nothing,
)
    isfile(database) || throw(ArgumentError("PHREEQC database not found: $database"))
    run_input = isfile(input) ? LibPhreeqcRM.RM_RunFile : LibPhreeqcRM.RM_RunString
    isfile(input) || occursin('\n', input) || throw(ArgumentError("PHREEQC input file not found: $input"))
    ic = _initial_conditions(initial, ncells)

    chemistry = Chemistry(LibPhreeqcRM.RM_Create(ncells, threads), ncells)
    try
        _initialize!(chemistry, database, input, run_input, ic, output)
    catch
        close(chemistry)
        rethrow()
    end
    return chemistry
end

function _initialize!(chemistry::Chemistry, database, input, run_input, ic, output)
    (; id, ncells) = chemistry
    rm = LibPhreeqcRM
    rm.RM_SetScreenOn(id, false)
    if output !== nothing
        rm.RM_SetFilePrefix(id, output)
        rm.RM_OpenFiles(id)
    end
    rm.RM_LoadDatabase(id, database)
    # Concentrations in mol/L; other reactants per litre of water; 1 L of water per cell
    rm.RM_SetUnitsSolution(id, 2)
    for set_units in (
        rm.RM_SetUnitsPPassemblage,
        rm.RM_SetUnitsExchange,
        rm.RM_SetUnitsSurface,
        rm.RM_SetUnitsGasPhase,
        rm.RM_SetUnitsSSassemblage,
        rm.RM_SetUnitsKinetics,
    )
        set_units(id, 1)
    end
    rm.RM_SetPorosity(id, ones(ncells))
    rm.RM_SetSaturationUser(id, ones(ncells))
    # Concentrations are moles per litre of water in the cell, which conserves the mass of
    # the transported species (the solution volume computed by PHREEQC is not used)
    rm.RM_UseSolutionDensityVolume(id, false)
    rm.RM_SetSpeciesSaveOn(id, true)
    rm.RM_SetSelectedOutputOn(id, true)

    # The input is also run in the worker instances, which compute the cells and need the
    # selected output and rate definitions; their solutions and reactants are then deleted,
    # so that the cells contain only the definitions assigned by `initial`.
    run_input(id, true, true, false, input)
    rm.RM_RunString(id, true, false, false, "DELETE; -all")
    ncomps = rm.RM_FindComponents(id)
    chemistry.components = [rm.RM_GetComponent(id, i) for i in 0:(ncomps - 1)]
    chemistry.species = [rm.RM_GetSpeciesName(id, i) for i in 0:(rm.RM_GetSpeciesCount(id) - 1)]
    rm.RM_InitialPhreeqc2Module(id, ic)
    rm.RM_SetTime(id, 0.0)
    rm.RM_SetTimeStep(id, 0.0)
    rm.RM_RunCells(id)
    return nothing
end

# (ncells × 7) matrix of definition numbers, -1 where an entity is absent
function _initial_conditions(initial::NamedTuple, ncells::Integer)
    for key in keys(initial)
        key in INITIAL_KEYS ||
            throw(ArgumentError("Unknown initial condition `$key`; the keys are $(join(INITIAL_KEYS, ", "))"))
    end
    ic = fill(-1, ncells, length(INITIAL_KEYS))
    for (j, key) in enumerate(INITIAL_KEYS)
        haskey(initial, key) || continue
        value = initial[key]
        value isa Integer ||
            length(value) == ncells ||
            throw(ArgumentError("`$key` must be a number or a vector of $ncells numbers"))
        ic[:, j] .= value
    end
    return ic
end

function Base.close(chemistry::Chemistry)
    if chemistry.id >= 0
        LibPhreeqcRM.RM_Destroy(chemistry.id)
        chemistry.id = -1
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", c::Chemistry)
    print(io, "Chemistry: ", c.ncells, " cells, ", length(c.components), " components, ")
    print(io, length(c.species), " aqueous species, t = ", c.time, " s")
    c.id < 0 && print(io, " (closed)")
    return nothing
end

"""
    aqueous_species(chemistry) -> Vector{Species}

Aqueous species of the PHREEQC database that are present in the system, with their charge
numbers and their diffusion coefficients at 25 °C (`-dw` in the database; zero if not
defined). The order is that of the columns of [`species_concentrations`](@ref).
"""
function aqueous_species(chemistry::Chemistry)
    z = LibPhreeqcRM.RM_GetSpeciesZ(chemistry.id)
    D = LibPhreeqcRM.RM_GetSpeciesD25(chemistry.id)
    return [Species(name, round(Int, z[i]), D[i]) for (i, name) in enumerate(chemistry.species)]
end

"""
    species_concentrations(chemistry) -> Matrix{Float64}

Concentrations of the aqueous species (mol m⁻³), as a (cells × species) matrix in the order
of [`aqueous_species`](@ref).
"""
species_concentrations(chemistry::Chemistry) = 1000.0 .* LibPhreeqcRM.RM_GetSpeciesConcentrations(chemistry.id)

"""
    set_species_concentrations!(chemistry, c)

Set the concentrations of the aqueous species (mol m⁻³) in all cells, given as a
(cells × species) matrix in the order of [`aqueous_species`](@ref). The solution of each
cell is redefined by the total amounts of the components.
"""
function set_species_concentrations!(chemistry::Chemistry, c::AbstractMatrix{<:Real})
    LibPhreeqcRM.RM_SpeciesConcentrations2Module(chemistry.id, c ./ 1000.0)
    return chemistry
end

"""
    component_concentrations(chemistry) -> Matrix{Float64}

Total concentrations of the components (mol m⁻³), as a (cells × components) matrix in the
order of `chemistry.components`.
"""
component_concentrations(chemistry::Chemistry) = 1000.0 .* LibPhreeqcRM.RM_GetConcentrations(chemistry.id)

"""
    react!(chemistry, Δt) -> chemistry

Compute the reactions in all cells over a time step `Δt` (s): kinetic reactions are
integrated over `Δt`, and equilibrium is then established. With `Δt = 0`, the cells are
only equilibrated.
"""
function react!(chemistry::Chemistry, Δt::Real)
    Δt >= 0 || throw(ArgumentError("The time step must be non-negative, got $Δt"))
    chemistry.time += Δt
    LibPhreeqcRM.RM_SetTime(chemistry.id, chemistry.time)
    LibPhreeqcRM.RM_SetTimeStep(chemistry.id, Δt)
    LibPhreeqcRM.RM_RunCells(chemistry.id)
    return chemistry
end

"""
    selected_output(chemistry; n_user=1) -> Dict{String, Vector{Float64}}

Columns of the SELECTED_OUTPUT block `n_user`, with one value per cell, indexed by their
headings (e.g. `"pH"`, `"si_Calcite"`).
"""
function selected_output(chemistry::Chemistry; n_user::Integer=1)
    id = chemistry.id
    LibPhreeqcRM.RM_SetCurrentSelectedOutputUserNumber(id, n_user)
    so = LibPhreeqcRM.RM_GetSelectedOutput(id)
    headings = [LibPhreeqcRM.RM_GetSelectedOutputHeading(id, j) for j in 0:(size(so, 2) - 1)]
    return Dict(heading => so[:, j] for (j, heading) in enumerate(headings))
end
