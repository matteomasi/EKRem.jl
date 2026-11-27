"""
    LibPhreeqcRM

Low-level interface to the C API of PhreeqcRM (`RM_interface_C.h`). Each C function is
represented by one Julia function with the same name and arguments. Functions that return
an `IRM_RESULT` status throw a [`PhreeqcRMError`](@ref) on failure, with the message
returned by `RM_GetErrorString`.

Only the functions required by SimEK are wrapped. The shared library is selected with
[`set_phreeqcrm_library`](@ref); see [`phreeqcrm_library`](@ref) for the order of precedence.
"""
module LibPhreeqcRM

using Libdl
using Preferences: load_preference, set_preferences!

const IRM_RESULT_NAMES = Dict{Cint,String}(
    0 => "IRM_OK",
    -1 => "IRM_OUTOFMEMORY",
    -2 => "IRM_BADVARTYPE",
    -3 => "IRM_INVALIDARG",
    -4 => "IRM_INVALIDROW",
    -5 => "IRM_INVALIDCOL",
    -6 => "IRM_BADINSTANCE",
    -7 => "IRM_FAIL",
)

"Error raised when a PhreeqcRM function returns a negative status."
struct PhreeqcRMError <: Exception
    func::Symbol
    status::Cint
    message::String
end

function Base.showerror(io::IO, e::PhreeqcRMError)
    name = get(IRM_RESULT_NAMES, e.status, string(e.status))
    print(io, "PhreeqcRM: ", e.func, " returned ", name)
    isempty(e.message) || print(io, "\n", rstrip(e.message))
    return nothing
end

# ------------------------------------------------------------------------------------------
# Shared library

const ENV_KEY = "SIMEK_PHREEQCRM"
const PREFERENCE_KEY = "phreeqcrm_library"

const _handle = Ref{Ptr{Cvoid}}(C_NULL)
const _path = Ref{String}("")

_package() = parentmodule(@__MODULE__)   # SimEK, owner of the preference

"""
    set_phreeqcrm_library(path; save=true)

Load the PhreeqcRM shared library (`PhreeqcRM.dll`, `libphreeqcrm.so`, ...) located at
`path` and return its absolute path.

With `save=true`, the path is also saved as the preference `phreeqcrm_library` in the file
`LocalPreferences.toml` of the active Julia project. The library is then found
automatically in later sessions and by `Pkg.test`, so the function needs to be called only
once per project (or once in the global environment). With `save=false`, the library is
used in the current session only.
"""
function set_phreeqcrm_library(path::AbstractString; save::Bool=true)
    file = abspath(path)
    _open(file)
    if save
        _save_library(file)
        @info "PhreeqcRM library saved as a preference" library = file project = Base.active_project()
    end
    return file
end

"""
    phreeqcrm_library() -> Union{String, Nothing}

Path of the PhreeqcRM library used by SimEK, or `nothing` if no library has been set. The
library is selected in the following order of precedence:

1. the library loaded in the current session by [`set_phreeqcrm_library`](@ref);
2. the environment variable `SIMEK_PHREEQCRM`;
3. the preference saved by `set_phreeqcrm_library` in a `LocalPreferences.toml` file of
   the active project or of another environment in the load path.
"""
function phreeqcrm_library()
    isempty(_path[]) || return _path[]
    path = get(ENV, ENV_KEY, "")
    isempty(path) || return abspath(path)
    return load_preference(_package(), PREFERENCE_KEY, nothing)
end

function _save_library(file::AbstractString; project_toml::AbstractString=Base.active_project())
    set_preferences!(_package(), PREFERENCE_KEY => file; force=true, project_toml=project_toml)
    return nothing
end

function _open(file::AbstractString)
    isfile(file) || error("PhreeqcRM library not found: $file")
    _handle[] = Libdl.dlopen(file)
    _path[] = file
    return nothing
end

function _library()
    if _handle[] == C_NULL
        file = phreeqcrm_library()
        file === nothing && error(
            "The PhreeqcRM library has not been set. Call `set_phreeqcrm_library(path)` once " *
            "(the path is saved for later sessions) or set the environment variable $ENV_KEY.",
        )
        _open(file)
    end
    return _handle[]
end

_fptr(name::Symbol) = Libdl.dlsym(_library(), name)

function _check(id::Integer, status::Integer, func::Symbol)
    status >= 0 && return nothing
    return throw(PhreeqcRMError(func, Cint(status), _error_string(id)))
end

function _error_string(id::Integer)
    l = @ccall $(_fptr(:RM_GetErrorStringLength))(id::Cint)::Cint
    l <= 0 && return ""
    buf = Vector{UInt8}(undef, l + 1)
    @ccall $(_fptr(:RM_GetErrorString))(id::Cint, buf::Ptr{UInt8}, (l + 1)::Cint)::Cint
    return unsafe_string(pointer(buf))
end

function _string_out(f::Function, len::Integer)
    buf = zeros(UInt8, len)
    f(buf)
    return unsafe_string(pointer(buf))
end

# ------------------------------------------------------------------------------------------
# Creation, configuration and files

"""
    RM_Create(nxyz, nthreads) -> id

Create a PhreeqcRM instance with `nxyz` cells and `nthreads` OpenMP threads
(`nthreads ≤ 0` selects the number of available processors). Returns the instance id.
"""
function RM_Create(nxyz::Integer, nthreads::Integer)
    id = @ccall $(_fptr(:RM_Create))(nxyz::Cint, nthreads::Cint)::Cint
    id < 0 && throw(PhreeqcRMError(:RM_Create, id, ""))
    return Int(id)
end

function RM_Destroy(id::Integer)
    s = @ccall $(_fptr(:RM_Destroy))(id::Cint)::Cint
    s < 0 && throw(PhreeqcRMError(:RM_Destroy, s, ""))  # the instance no longer exists
    return nothing
end

function RM_LoadDatabase(id::Integer, db_name::AbstractString)
    s = @ccall $(_fptr(:RM_LoadDatabase))(id::Cint, db_name::Cstring)::Cint
    return _check(id, s, :RM_LoadDatabase)
end

function RM_SetFilePrefix(id::Integer, prefix::AbstractString)
    s = @ccall $(_fptr(:RM_SetFilePrefix))(id::Cint, prefix::Cstring)::Cint
    return _check(id, s, :RM_SetFilePrefix)
end

function RM_OpenFiles(id::Integer)
    s = @ccall $(_fptr(:RM_OpenFiles))(id::Cint)::Cint
    return _check(id, s, :RM_OpenFiles)
end

function RM_CloseFiles(id::Integer)
    s = @ccall $(_fptr(:RM_CloseFiles))(id::Cint)::Cint
    return _check(id, s, :RM_CloseFiles)
end

"Print (or not) PhreeqcRM messages on the screen; errors are always available from `RM_GetErrorString`."
function RM_SetScreenOn(id::Integer, tf::Bool)
    s = @ccall $(_fptr(:RM_SetScreenOn))(id::Cint, tf::Cint)::Cint
    return _check(id, s, :RM_SetScreenOn)
end

"Include (or not) H2O as a component, in addition to the total H and O."
function RM_SetComponentH2O(id::Integer, tf::Bool)
    s = @ccall $(_fptr(:RM_SetComponentH2O))(id::Cint, tf::Cint)::Cint
    return _check(id, s, :RM_SetComponentH2O)
end

"""
    RM_UseSolutionDensityVolume(id, tf)

Convert between moles and concentrations with the solution volume computed by PHREEQC
(`true`, default of PhreeqcRM) or with the volume porosity × saturation × representative
volume of the cell (`false`).
"""
function RM_UseSolutionDensityVolume(id::Integer, tf::Bool)
    s = @ccall $(_fptr(:RM_UseSolutionDensityVolume))(id::Cint, tf::Cint)::Cint
    return _check(id, s, :RM_UseSolutionDensityVolume)
end

function _set_option(id::Integer, func::Symbol, option::Integer)
    s = @ccall $(_fptr(func))(id::Cint, option::Cint)::Cint
    return _check(id, s, func)
end

# Units of the solution: 1 = mg/L, 2 = mol/L, 3 = kg/kgs. Units of the other reactants:
# 0 = mol/L of cell, 1 = mol/L of water, 2 = mol/L of rock.
RM_SetUnitsSolution(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsSolution, option)
RM_SetUnitsPPassemblage(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsPPassemblage, option)
RM_SetUnitsExchange(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsExchange, option)
RM_SetUnitsSurface(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsSurface, option)
RM_SetUnitsGasPhase(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsGasPhase, option)
RM_SetUnitsSSassemblage(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsSSassemblage, option)
RM_SetUnitsKinetics(id::Integer, option::Integer) = _set_option(id, :RM_SetUnitsKinetics, option)

function _set_cell_values(id::Integer, func::Symbol, values::AbstractVector{<:Real})
    length(values) == RM_GetGridCellCount(id) || throw(ArgumentError("$func: one value per cell is required"))
    v = Vector{Cdouble}(values)
    s = @ccall $(_fptr(func))(id::Cint, v::Ptr{Cdouble})::Cint
    return _check(id, s, func)
end

# Properties of the cells, one value per cell: porosity, saturation, and representative
# volume (L)
RM_SetPorosity(id::Integer, por::AbstractVector{<:Real}) = _set_cell_values(id, :RM_SetPorosity, por)
RM_SetSaturationUser(id::Integer, sat::AbstractVector{<:Real}) = _set_cell_values(id, :RM_SetSaturationUser, sat)
function RM_SetRepresentativeVolume(id::Integer, rv::AbstractVector{<:Real})
    return _set_cell_values(id, :RM_SetRepresentativeVolume, rv)
end

# ------------------------------------------------------------------------------------------
# Running PHREEQC input

"""
    RM_RunFile(id, workers, initial_phreeqc, utility, chem_name)

Run the PHREEQC input file `chem_name` in the worker, InitialPhreeqc and/or Utility
instances, as selected by the three Boolean flags.
"""
function RM_RunFile(id::Integer, workers::Bool, initial_phreeqc::Bool, utility::Bool, chem_name::AbstractString)
    s = @ccall $(_fptr(:RM_RunFile))(
        id::Cint, workers::Cint, initial_phreeqc::Cint, utility::Cint, chem_name::Cstring
    )::Cint
    return _check(id, s, :RM_RunFile)
end

"""
    RM_RunString(id, workers, initial_phreeqc, utility, input_string)

Run PHREEQC input given as a string, in the instances selected by the three Boolean flags.
"""
function RM_RunString(id::Integer, workers::Bool, initial_phreeqc::Bool, utility::Bool, input_string::AbstractString)
    s = @ccall $(_fptr(:RM_RunString))(
        id::Cint, workers::Cint, initial_phreeqc::Cint, utility::Cint, input_string::Cstring
    )::Cint
    return _check(id, s, :RM_RunString)
end

# ------------------------------------------------------------------------------------------
# Components and grid

"""
    RM_FindComponents(id) -> ncomps

Determine the components (elements, charge balance, and optionally H, O, H2O) from all
definitions run so far. Returns the number of components.
"""
function RM_FindComponents(id::Integer)
    n = @ccall $(_fptr(:RM_FindComponents))(id::Cint)::Cint
    _check(id, n, :RM_FindComponents)
    return Int(n)
end

"""
    RM_GetComponent(id, num) -> name

Name of component `num` (zero-based, as in the C interface).
"""
function RM_GetComponent(id::Integer, num::Integer; maxlen::Integer=256)
    return _string_out(maxlen) do buf
        s = @ccall $(_fptr(:RM_GetComponent))(id::Cint, num::Cint, buf::Ptr{UInt8}, maxlen::Cint)::Cint
        return _check(id, s, :RM_GetComponent)
    end
end

"Number of components found by [`RM_FindComponents`](@ref)."
function RM_GetComponentCount(id::Integer)
    n = @ccall $(_fptr(:RM_GetComponentCount))(id::Cint)::Cint
    _check(id, n, :RM_GetComponentCount)
    return Int(n)
end

function RM_GetGridCellCount(id::Integer)
    n = @ccall $(_fptr(:RM_GetGridCellCount))(id::Cint)::Cint
    _check(id, n, :RM_GetGridCellCount)
    return Int(n)
end

# ------------------------------------------------------------------------------------------
# Initial conditions and reaction step

"""
    RM_InitialPhreeqc2Module(id, ic1)

Transfer solutions and reactants defined in the InitialPhreeqc instance to the cells of the
module. `ic1` is an (nxyz, 7) integer matrix of definition numbers, with the columns
solution, equilibrium phases, exchange, surface, gas phase, solid solutions and kinetics;
negative values mean that the entity is not defined in the cell. Mixing with a second set of
definitions is not used (`initial_conditions2` and `fraction1` are NULL in the C call).
"""
function RM_InitialPhreeqc2Module(id::Integer, ic1::AbstractMatrix{<:Integer})
    nxyz = RM_GetGridCellCount(id)
    size(ic1) == (nxyz, 7) || throw(ArgumentError("ic1 must be a $nxyz × 7 matrix, got $(size(ic1))"))
    c = Matrix{Cint}(ic1)
    s = @ccall $(_fptr(:RM_InitialPhreeqc2Module))(
        id::Cint, c::Ptr{Cint}, C_NULL::Ptr{Cint}, C_NULL::Ptr{Cdouble}
    )::Cint
    return _check(id, s, :RM_InitialPhreeqc2Module)
end

"Run the reactions (equilibrium and kinetics over the current time step) in all cells."
function RM_RunCells(id::Integer)
    s = @ccall $(_fptr(:RM_RunCells))(id::Cint)::Cint
    return _check(id, s, :RM_RunCells)
end

"Set the simulation time (s) at the end of the reaction step."
function RM_SetTime(id::Integer, time::Real)
    s = @ccall $(_fptr(:RM_SetTime))(id::Cint, time::Cdouble)::Cint
    return _check(id, s, :RM_SetTime)
end

"Set the length of the reaction step (s), over which the kinetic reactions are integrated."
function RM_SetTimeStep(id::Integer, time_step::Real)
    s = @ccall $(_fptr(:RM_SetTimeStep))(id::Cint, time_step::Cdouble)::Cint
    return _check(id, s, :RM_SetTimeStep)
end

"""
    RM_GetConcentrations(id) -> Matrix{Float64}

Concentrations of all components in all cells, as a (cells × components) matrix, in the
units set by `RM_SetUnitsSolution`. The array passed to the C function is allocated here.
"""
function RM_GetConcentrations(id::Integer)
    c = zeros(Float64, RM_GetGridCellCount(id), RM_GetComponentCount(id))
    s = @ccall $(_fptr(:RM_GetConcentrations))(id::Cint, c::Ptr{Cdouble})::Cint
    _check(id, s, :RM_GetConcentrations)
    return c
end

# ------------------------------------------------------------------------------------------
# Selected output

"Enable or disable the selected output of the worker instances (required by `RM_GetSelectedOutput`)."
function RM_SetSelectedOutputOn(id::Integer, selected_output::Bool)
    s = @ccall $(_fptr(:RM_SetSelectedOutputOn))(id::Cint, selected_output::Cint)::Cint
    return _check(id, s, :RM_SetSelectedOutputOn)
end

"Number of SELECTED_OUTPUT blocks defined in the worker instances."
function RM_GetSelectedOutputCount(id::Integer)
    n = @ccall $(_fptr(:RM_GetSelectedOutputCount))(id::Cint)::Cint
    _check(id, n, :RM_GetSelectedOutputCount)
    return Int(n)
end

"Select the SELECTED_OUTPUT block (user number `n_user`) returned by the other selected-output functions."
function RM_SetCurrentSelectedOutputUserNumber(id::Integer, n_user::Integer)
    s = @ccall $(_fptr(:RM_SetCurrentSelectedOutputUserNumber))(id::Cint, n_user::Cint)::Cint
    return _check(id, s, :RM_SetCurrentSelectedOutputUserNumber)
end

function RM_GetSelectedOutputColumnCount(id::Integer)
    n = @ccall $(_fptr(:RM_GetSelectedOutputColumnCount))(id::Cint)::Cint
    _check(id, n, :RM_GetSelectedOutputColumnCount)
    return Int(n)
end

function RM_GetSelectedOutputRowCount(id::Integer)
    n = @ccall $(_fptr(:RM_GetSelectedOutputRowCount))(id::Cint)::Cint
    _check(id, n, :RM_GetSelectedOutputRowCount)
    return Int(n)
end

"""
    RM_GetSelectedOutputHeading(id, icol) -> heading

Heading of column `icol` (zero-based, as in the C interface) of the current selected output.
"""
function RM_GetSelectedOutputHeading(id::Integer, icol::Integer; maxlen::Integer=256)
    return _string_out(maxlen) do buf
        s = @ccall $(_fptr(:RM_GetSelectedOutputHeading))(id::Cint, icol::Cint, buf::Ptr{UInt8}, maxlen::Cint)::Cint
        return _check(id, s, :RM_GetSelectedOutputHeading)
    end
end

"""
    RM_GetSelectedOutput(id) -> Matrix{Float64}

Current selected output as a (cells × columns) matrix. The array passed to the C function
is allocated here, with the cell index varying fastest as in PhreeqcRM.
"""
function RM_GetSelectedOutput(id::Integer)
    so = zeros(Float64, RM_GetSelectedOutputRowCount(id), RM_GetSelectedOutputColumnCount(id))
    s = @ccall $(_fptr(:RM_GetSelectedOutput))(id::Cint, so::Ptr{Cdouble})::Cint
    _check(id, s, :RM_GetSelectedOutput)
    return so
end

# ------------------------------------------------------------------------------------------
# Aqueous species (transport of individual species)

"""
    RM_SetSpeciesSaveOn(id, save_on)

Enable or disable the storage of aqueous-species information, which is required for the
transport of individual species. Must be called before [`RM_FindComponents`](@ref).
"""
function RM_SetSpeciesSaveOn(id::Integer, save_on::Bool)
    s = @ccall $(_fptr(:RM_SetSpeciesSaveOn))(id::Cint, save_on::Cint)::Cint
    return _check(id, s, :RM_SetSpeciesSaveOn)
end

"Number of aqueous species found by [`RM_FindComponents`](@ref) (species saving enabled)."
function RM_GetSpeciesCount(id::Integer)
    n = @ccall $(_fptr(:RM_GetSpeciesCount))(id::Cint)::Cint
    _check(id, n, :RM_GetSpeciesCount)
    return Int(n)
end

"""
    RM_GetSpeciesName(id, i) -> name

Name of aqueous species `i` (zero-based, as in the C interface).
"""
function RM_GetSpeciesName(id::Integer, i::Integer; maxlen::Integer=256)
    return _string_out(maxlen) do buf
        s = @ccall $(_fptr(:RM_GetSpeciesName))(id::Cint, i::Cint, buf::Ptr{UInt8}, maxlen::Cint)::Cint
        return _check(id, s, :RM_GetSpeciesName)
    end
end

"""
    RM_GetSpeciesZ(id) -> Vector{Float64}

Charge numbers of all aqueous species, in the order of [`RM_GetSpeciesName`](@ref). The
vector passed to the C function is allocated here.
"""
function RM_GetSpeciesZ(id::Integer)
    z = zeros(Float64, RM_GetSpeciesCount(id))
    s = @ccall $(_fptr(:RM_GetSpeciesZ))(id::Cint, z::Ptr{Cdouble})::Cint
    _check(id, s, :RM_GetSpeciesZ)
    return z
end

"""
    RM_GetSpeciesConcentrations(id) -> Matrix{Float64}

Concentrations (mol L⁻¹) of all aqueous species in all cells, as a (cells × species) matrix
in the order of [`RM_GetSpeciesName`](@ref). Inactive cells are set to 1e30 by PhreeqcRM.
The array passed to the C function is allocated here.
"""
function RM_GetSpeciesConcentrations(id::Integer)
    c = zeros(Float64, RM_GetGridCellCount(id), RM_GetSpeciesCount(id))
    s = @ccall $(_fptr(:RM_GetSpeciesConcentrations))(id::Cint, c::Ptr{Cdouble})::Cint
    _check(id, s, :RM_GetSpeciesConcentrations)
    return c
end

"""
    RM_SpeciesConcentrations2Module(id, c)

Transfer the concentrations (mol L⁻¹) of all aqueous species in all cells, given as a
(cells × species) matrix, to the module. The component totals of each cell are computed
from the species concentrations.
"""
function RM_SpeciesConcentrations2Module(id::Integer, c::AbstractMatrix{<:Real})
    dims = (RM_GetGridCellCount(id), RM_GetSpeciesCount(id))
    size(c) == dims || throw(ArgumentError("The concentrations must be a $(dims[1]) × $(dims[2]) matrix"))
    m = Matrix{Cdouble}(c)
    s = @ccall $(_fptr(:RM_SpeciesConcentrations2Module))(id::Cint, m::Ptr{Cdouble})::Cint
    return _check(id, s, :RM_SpeciesConcentrations2Module)
end

"""
    RM_GetSolutionVolume(id) -> Vector{Float64}

Volume of the solution (L) in each cell, as computed by PHREEQC. The vector passed to the C
function is allocated here.
"""
function RM_GetSolutionVolume(id::Integer)
    v = zeros(Float64, RM_GetGridCellCount(id))
    s = @ccall $(_fptr(:RM_GetSolutionVolume))(id::Cint, v::Ptr{Cdouble})::Cint
    _check(id, s, :RM_GetSolutionVolume)
    return v
end

"""
    RM_GetSpeciesD25(id) -> Vector{Float64}

Diffusion coefficients at 25 °C (m² s⁻¹) of all aqueous species, as defined by `-dw` in
the database, in the order of [`RM_GetSpeciesName`](@ref). Species without `-dw` have a
coefficient of zero. The vector passed to the C function is allocated here.
"""
function RM_GetSpeciesD25(id::Integer)
    d = zeros(Float64, RM_GetSpeciesCount(id))
    s = @ccall $(_fptr(:RM_GetSpeciesD25))(id::Cint, d::Ptr{Cdouble})::Cint
    _check(id, s, :RM_GetSpeciesD25)
    return d
end

end # module LibPhreeqcRM
