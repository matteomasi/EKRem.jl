# One-dimensional model: columns and laboratory cells.

@doc raw"""
    Model1D(; length, grid, species, initial, left, right, potential=Poisson(),
            porosity=1.0, tortuosity=1.0, relative_permittivity=WATER_RELATIVE_PERMITTIVITY,
            electroosmotic_permeability=0.0, temperature=298.15, area=nothing,
            database=nothing, input=nothing, threads=1)

One-dimensional model of electrokinetic transport in a water-saturated porous medium
between two boundaries at ``x = 0`` (`left`) and ``x = L`` (`right`), with ``L`` = `length`
(m).

# Arguments
- `grid`: node coordinates (a vector from 0 to `length`), or a named tuple
  `(h_min=…, h_max=…, growth=…)` for a grid refined toward both ends
  (see [`grid_refined`](@ref)).
- `species`: vector of [`Species`](@ref).
- `initial`: uniform initial pore-water concentrations (mol m⁻³), as a dictionary that maps
  species names to concentrations.
- `left`, `right`: boundary conditions, [`Reservoir`](@ref) or [`Compartment`](@ref). The
  potential must be fixed at one boundary at least; a current density can be applied at
  the other one with the [`ChargeConservation`](@ref) formulation.
- `potential`: formulation of the potential equation: [`Poisson`](@ref),
  [`ChargeConservation`](@ref) or [`FixedField`](@ref).
- `porosity`, `tortuosity`, `relative_permittivity`, `electroosmotic_permeability`:
  properties of the medium (see [`PorousMedium`](@ref)).
- `temperature`: temperature (K).
- `area`: cross-sectional area of the column (m²), required with electrode compartments.

# Electroosmosis
The electroosmotic water flux is uniform along the column, as required by the
conservation of the water volume. It is the mean of ``-k_{eo} \partial\phi/\partial x``
over the column (Masi et al., 2017, Eq. 8):

```math
q = \frac{k_{eo}}{L} \left[ \phi(0) - \phi(L) \right],
```

which corresponds to equal hydraulic heads at both ends. Water flows in from one boundary
and out at the other; a compartment downstream overflows, and a compartment upstream is
refilled with deionised water, so that their volumes remain constant.

The initial and boundary compositions must be electroneutral. The model is run with
[`run!`](@ref), and the results are obtained with [`concentration`](@ref),
[`electric_potential`](@ref), [`voltage`](@ref), [`current_density`](@ref) and
[`electroosmotic_flow`](@ref) at the times in [`output_times`](@ref).

# Geochemical reactions
Reactions are computed by PhreeqcRM in every node (see [`Chemistry`](@ref)) when the PHREEQC
`database` and `input` are given. Then:
- `initial` contains the numbers of the PHREEQC definitions in the column, e.g.
  `(solution=1, kinetics=1)`, as for [`Chemistry`](@ref), with one number for all nodes or
  one per node. The amounts of reactants are per litre of pore water.
- The reservoirs and compartments are given by PHREEQC definitions, e.g.
  `Reservoir(; solution=2, potential=0.0)`.
- `species` is optional: by default, all aqueous species except H2O are transported, with
  the diffusion coefficients of the database. A vector of names or of [`Species`](@ref)
  selects the transported species (and their diffusion coefficients); the other species
  are immobile.
- `threads` is the number of threads of PhreeqcRM.
"""
mutable struct Model1D
    x::Vector{Float64}
    species::Vector{Species}
    medium::PorousMedium
    temperature::Float64
    potential::PotentialModel
    left::BoundaryCondition
    right::BoundaryCondition
    problem::TransportProblem
    chemistry::Union{Coupling,Nothing}
    state::Matrix{Float64}
    time::Float64
    times::Vector{Float64}
    results::Vector{Matrix{Float64}}
end

function Model1D(;
    length::Real,
    grid,
    species=nothing,
    initial::Union{AbstractDict,NamedTuple},
    left::BoundaryCondition,
    right::BoundaryCondition,
    potential::PotentialModel=Poisson(),
    porosity::Real=1.0,
    tortuosity::Real=1.0,
    relative_permittivity::Real=WATER_RELATIVE_PERMITTIVITY,
    electroosmotic_permeability::Real=0.0,
    temperature::Real=298.15,
    area::Union{Real,Nothing}=nothing,
    database::Union{AbstractString,Nothing}=nothing,
    input::Union{AbstractString,Nothing}=nothing,
    threads::Integer=1,
)
    x = _coordinates(length, grid)
    medium = PorousMedium(; porosity, tortuosity, relative_permittivity, electroosmotic_permeability)
    _check_electric_conditions(potential, left, right)

    # Concentrations at the nodes (species × nodes), the boundary nodes being the reservoirs
    if database === nothing && input === nothing
        initial isa AbstractDict ||
            throw(ArgumentError("`initial` must be a dictionary of concentrations (no PHREEQC database and input)"))
        species isa AbstractVector{Species} || throw(ArgumentError("`species` must be a vector of Species"))
        any(b -> _phreeqc_boundary(b), (left, right)) &&
            throw(ArgumentError("Boundaries given by PHREEQC definitions require a PHREEQC `database` and `input`"))
        coupling = nothing
        c = repeat(concentration_vector(species, initial), 1, Base.length(x))
        c[:, 1] = concentration_vector(species, left.concentrations)
        c[:, end] = concentration_vector(species, right.concentrations)
    else
        coupling, species, c = _chemistry(x, species, initial, left, right, database, input, threads)
        left, right = (_resolved(b, species, c[:, k]) for (b, k) in ((left, 1), (right, Base.length(x))))
    end
    names = [s.name for s in species]
    allunique(names) || throw(ArgumentError("Species names must be unique"))
    for k in axes(c, 2)
        imbalance = charge_imbalance(species, c[:, k])
        label = if k == 1
            "left boundary"
        elseif k == Base.length(x)
            "right boundary"
        else
            "initial"
        end
        imbalance <= 1.0e-8 ||
            throw(ArgumentError("The $label composition is not electroneutral (relative charge imbalance $imbalance)"))
    end

    # Characteristic scales: largest concentration of each species (with a floor for species
    # absent from all compositions) and R T / F for the potential
    cmax = vec(maximum(c; dims=2))
    floor_c = 1.0e-9 * max(maximum(cmax), 1.0)
    scales = vcat(max.(cmax, floor_c), GAS_CONSTANT * temperature / FARADAY)

    electric(b) = fixed_potential(b) ? b.potential : b.current_density
    bvalues = hcat(vcat(c[:, 1], electric(left)), vcat(c[:, end], electric(right)))
    depth, flushing, production, electrode_current = _compartments(species, left, right, area)
    problem = TransportProblem(
        simplexgrid(x),
        species,
        medium,
        potential,
        bvalues,
        [fixed_potential(left), fixed_potential(right)],
        scales;
        temperature,
        depth,
        flushing,
        production,
        electrode_current,
    )

    # Initial state. The potential is computed from the potential equation, starting from a
    # linear profile (with an applied current density, the field is estimated from Ohm's law
    # with the mean initial conductivity).
    c_mean = vec(sum(c; dims=2)) ./ Base.length(x)
    φ0 = _initial_potential(x, left, right, problem.data.Deff, [s.z for s in species], c_mean, temperature)
    state = solve_potential(problem, vcat(c, transpose(φ0)))

    return Model1D(
        x,
        collect(species),
        medium,
        temperature,
        potential,
        left,
        right,
        problem,
        coupling,
        state,
        0.0,
        [0.0],
        [copy(state)],
    )
end

# PhreeqcRM cells at the nodes, transported species, and their concentrations at the nodes
function _chemistry(x, species, initial, left, right, database, input, threads)
    database !== nothing && input !== nothing ||
        throw(ArgumentError("Geochemical reactions require both a PHREEQC `database` and an `input`"))
    initial isa NamedTuple ||
        throw(ArgumentError("With geochemical reactions, `initial` contains the numbers of PHREEQC definitions"))
    nnodes = Base.length(x)
    conditions = node_conditions(initial, nnodes, phreeqc_definitions(left), phreeqc_definitions(right))
    chemistry = Chemistry(nnodes; database, input, initial=conditions, threads)
    transported, columns = transported_species(chemistry, species)
    c = Matrix(transpose(species_concentrations(chemistry)[:, columns]))
    # Nodes whose concentrations change by the reactions: the interior and the compartments
    nodes = [k for k in 1:nnodes if !((k == 1 && left isa Reservoir) || (k == nnodes && right isa Reservoir))]
    coupling = Coupling(chemistry, columns, nodes, [_selected_output(chemistry)])
    return coupling, transported, c
end

_phreeqc_boundary(b::Reservoir) = b.solution !== nothing
_phreeqc_boundary(b::Compartment) = b.initial !== nothing

# Boundary condition with the initial concentrations computed by PhreeqcRM
function _resolved(b::Reservoir, species, c)
    return Reservoir(Dict(s.name => ci for (s, ci) in zip(species, c)), b.solution, b.potential, b.current_density)
end

function _resolved(b::Compartment, species, c)
    concentrations = Dict(s.name => ci for (s, ci) in zip(species, c))
    return Compartment(concentrations, b.initial, b.volume, b.flow, b.potential, b.current_density)
end

# Compartment volume and flushing rate per cross-sectional area, and production by water
# electrolysis (species × boundary). Under current control, the current density i applied at
# one boundary is positive where it flows into the column (anode); the current is the same at
# the other boundary (cathode). Under voltage control, the production follows from the
# computed current (`electrode_current`), and both H+ and OH- are required.
function _compartments(species, left, right, area)
    depth, flushing, production = zeros(2), zeros(2), zeros(Base.length(species), 2)
    electrode_current = falses(2)
    boundaries = (left, right)
    any(b -> b isa Compartment, boundaries) || return depth, flushing, production, electrode_current
    area isa Real && area > 0 ||
        throw(ArgumentError("Electrode compartments require the cross-sectional `area` of the column (m²)"))
    function index(product)
        j = findfirst(s -> s.name == product, species)
        j === nothing && throw(ArgumentError("Water electrolysis in a compartment requires the species $product"))
        return j
    end
    applied = findfirst(b -> !fixed_potential(b), boundaries)
    for (k, b) in enumerate(boundaries)
        b isa Compartment || continue
        depth[k], flushing[k] = b.volume / area, b.flow / area
        if applied === nothing
            index("H+"), index("OH-")
            electrode_current[k] = true
        else
            i = boundaries[applied].current_density
            anode = i >= 0 ? applied : 3 - applied
            production[index(k == anode ? "H+" : "OH-"), k] = abs(i) / FARADAY
        end
    end
    return depth, flushing, production, electrode_current
end

function _check_electric_conditions(potential::PotentialModel, left::BoundaryCondition, right::BoundaryCondition)
    fixed_potential(left) ||
        fixed_potential(right) ||
        throw(ArgumentError("The potential must be fixed at one boundary at least"))
    fixed_potential(left) && fixed_potential(right) && return nothing
    potential isa ChargeConservation || throw(
        ArgumentError(
            "An applied current density requires the ChargeConservation formulation, not $(nameof(typeof(potential)))",
        ),
    )
    return nothing
end

function _initial_potential(x, left::BoundaryCondition, right::BoundaryCondition, Deff, z, c0, temperature)
    L = x[end]
    fixed_potential(left) &&
        fixed_potential(right) &&
        return left.potential .+ (right.potential - left.potential) .* x ./ L
    σ = conductivity(c0, z, Deff; T=temperature)
    σ > 0 || throw(ArgumentError("An applied current density requires a conducting initial solution"))
    # A current density applied at the left boundary flows in the +x direction, and one
    # applied at the right boundary in the −x direction
    fixed_potential(right) && return right.potential .+ left.current_density / σ .* (L .- x)
    return left.potential .+ right.current_density / σ .* x
end

function _coordinates(length::Real, grid::NamedTuple)
    return grid_refined(length; grid...)
end

function _coordinates(length::Real, grid::AbstractVector)
    x = collect(Float64, grid)
    issorted(x) && allunique(x) || throw(ArgumentError("Grid coordinates must be strictly increasing"))
    tol = 1.0e-9 * length
    (abs(x[1]) <= tol && abs(x[end] - length) <= tol) ||
        throw(ArgumentError("Grid coordinates must range from 0 to $length"))
    x[1], x[end] = 0.0, Float64(length)
    return x
end

"""
    debye_ratio(model) -> Float64

Largest Debye length of the initial and boundary compositions divided by the smallest grid
spacing. For values much smaller than 1, the Debye layers are not resolved by the grid and
the solution of the Poisson equation is electroneutral outside thin boundary layers.
"""
function debye_ratio(model::Model1D)
    z = [s.z for s in model.species]
    n = Base.length(model.species)
    c = model.results[1][1:n, :]
    λ = maximum(
        debye_length(c[:, k], z; εr=model.medium.relative_permittivity, T=model.temperature) for
        k in axes(c, 2) if any(>(0), abs.(z) .* c[:, k])
    )
    return λ / minimum(diff(model.x))
end

function Base.show(io::IO, ::MIME"text/plain", m::Model1D)
    println(io, "Model1D: ", Base.length(m.x), " nodes on [0, ", m.x[end], "] m, h_min = ", minimum(diff(m.x)), " m")
    println(io, "  species: ", join((s.name for s in m.species), ", "))
    println(
        io,
        "  porosity = ",
        m.medium.porosity,
        ", tortuosity = ",
        m.medium.tortuosity,
        ", εr = ",
        m.medium.relative_permittivity,
        ", T = ",
        m.temperature,
        " K",
    )
    function electric(b)
        condition = fixed_potential(b) ? "φ = $(b.potential) V" : "i = $(b.current_density) A/m²"
        return b isa Compartment ? "compartment of $(b.volume) m³, $condition" : condition
    end
    println(io, "  boundaries: left ", electric(m.left), ", right ", electric(m.right))
    k_eo = m.medium.electroosmotic_permeability
    k_eo == 0 || println(io, "  electroosmosis: k_eo = ", k_eo, " m²/(V s)")
    m.chemistry === nothing ||
        println(io, "  reactions: PhreeqcRM, ", Base.length(m.chemistry.chemistry.components), " components")
    return print(
        io, "  potential: ", nameof(typeof(m.potential)), "; λ_D,max / h_min = ", round(debye_ratio(m); sigdigits=3)
    )
end

"""
    run!(model, times; coupling_step=nothing, solver=:Rodas5P, reltol=1e-4, abstol=1e-6,
         progress=true) -> model

Run the model from its current time to the largest of `times` (s), storing the
solution at each point of `times` vector. Repeated calls continue from the last state.

The time step is adapted automatically. `reltol` is the relative error tolerance of the
time integration of the concentrations, and `abstol` the absolute tolerance relative to the
largest concentration of each species. 

`progress=true` shows the simulated percentage and the elapsed time during the run. The
display is updated every 0.5 s by default; the interval is changed with, e.g.,
`EKRem.PROGRESS_INTERVAL[] = 5.0`.

`solver` selects the time integration method (both from OrdinaryDiffEq.jl, with error
control):
- `:Rodas5P` (default): Rosenbrock method of order 5. Accurate and fast for sharp fronts
  that move by electromigration; one Jacobian factorisation per step.
- `:FBDF`: backward differentiation formulas of variable order (1 to 5). Fewer Jacobian
  evaluations per step; suited to long simulations with slowly varying solutions.

Rodas5P conserves the amounts of the species to round-off; FBDF conserves them only to the
tolerance of its Newton iteration (relative errors of the order of `reltol`).

Solvers can be mixed for different parts of the simulation by calling `run!` multiple 
times with different `solver` arguments, e.g., see the example below.

With geochemical reactions, `coupling_step` (s) is required: each interval between output
times is divided into equal coupling steps not longer than `coupling_step`. Over each step,
the transport equations are integrated first, and the reactions are then computed over the
same step in every node (sequential non-iterative approach). The splitting error decreases
in proportion to the coupling step, and it is largest at sharp reaction fronts, such as an
acid front in a carbonate medium; the coupling step is therefore chosen by comparing the
results with those of a shorter step.

# Example
```julia
run!(model, [3600.0, 86400.0])                  # Rodas5P, output after 1 h and 1 d
run!(model, [7 * 86400.0]; solver=:FBDF)        # continue for 7 d with FBDF
```
"""
function run!(
    model::Model1D,
    times;
    coupling_step::Union{Real,Nothing}=nothing,
    solver::Symbol=:Rodas5P,
    reltol::Real=1.0e-4,
    abstol::Real=1.0e-6,
    progress::Bool=true,
)
    t = sort(unique(collect(Float64, times)))
    t[1] > model.time || throw(ArgumentError("Output times must be later than the current time $(model.time) s"))
    if model.chemistry === nothing
        coupling_step === nothing || throw(ArgumentError("`coupling_step` requires geochemical reactions"))
        states = solve_transport(model.problem, model.state, vcat(model.time, t); solver, reltol, abstol, progress)
        for (k, tk) in enumerate(t)
            push!(model.times, tk)
            push!(model.results, copy(states[k + 1]))
        end
    else
        coupling_step === nothing && throw(ArgumentError("Geochemical reactions require a `coupling_step`"))
        coupling_step > 0 || throw(ArgumentError("The coupling step must be positive"))
        _run_coupled!(model, t, coupling_step; solver, reltol, abstol, progress)
    end
    model.state = copy(model.results[end])
    model.time = model.times[end]
    return model
end

function _run_coupled!(model::Model1D, times, coupling_step; solver, reltol, abstol, progress)
    report = Progress(model.time, times[end])
    state = copy(model.state)
    t0 = model.time
    for tk in times
        nsteps = ceil(Int, (tk - t0) / coupling_step * (1 - 1.0e-12))
        steps = range(t0, tk; length=nsteps + 1)
        for i in 1:nsteps
            state = solve_transport(model.problem, state, steps[i:(i + 1)]; solver, reltol, abstol)[end]
            react_nodes!(model.chemistry, state, steps[i + 1] - steps[i])
            progress && show_progress!(report, steps[i + 1])
        end
        # The potential is computed for the concentrations after the reactions
        push!(model.times, tk)
        push!(model.results, solve_potential(model.problem, state, tk))
        push!(model.chemistry.selected, _selected_output(model.chemistry.chemistry))
        t0 = tk
    end
    progress && show_progress!(report, times[end]; final=true)
    return nothing
end

"Node coordinates (m) of the model."
coordinates(model::Model1D) = model.x

"Times (s) at which results are stored, starting with the initial state at t = 0."
output_times(model::Model1D) = model.times

function _time_index(model::Model1D, t::Real)
    k = findfirst(s -> isapprox(s, t; rtol=1.0e-12, atol=1.0e-9), model.times)
    k === nothing && throw(ArgumentError("No results at t = $t s; the output times are $(model.times)"))
    return k
end

"""
    concentration(model, name, t) -> Vector{Float64}

Pore-water concentration (mol m⁻³) of species `name` at the nodes, at output time `t` (s).
"""
function concentration(model::Model1D, name::AbstractString, t::Real)
    i = findfirst(s -> s.name == name, model.species)
    i === nothing && throw(ArgumentError("Unknown species \"$name\""))
    return model.results[_time_index(model, t)][i, :]
end

"""
    amount(model, name, t) -> Float64

Amount of species `name` (mol per m² of cross-sectional area) dissolved in the pore water of
the column and in the electrode compartments, at output time `t` (s).
"""
function amount(model::Model1D, name::AbstractString, t::Real)
    c = concentration(model, name, t)
    x = model.x
    volume = vcat(x[2] - x[1], x[3:end] .- x[1:(end - 2)], x[end] - x[end - 1]) ./ 2   # Voronoi cells
    d = model.problem.data
    compartments = (d.compartment[1] ? d.depth[1] * c[1] : 0.0) + (d.compartment[2] ? d.depth[2] * c[end] : 0.0)
    return model.medium.porosity * sum(volume .* c) + compartments
end

"""
    electric_potential(model, t) -> Vector{Float64}

Electric potential (V) at the nodes, at output time `t` (s).
"""
electric_potential(model::Model1D, t::Real) = model.results[_time_index(model, t)][end, :]

"""
    selected_output(model, t) -> Dict{String, Vector{Float64}}

Columns of the PHREEQC selected output at the nodes, indexed by their headings (e.g.
`"pH"`), at output time `t` (s). Only for models with geochemical reactions.
"""
function selected_output(model::Model1D, t::Real)
    model.chemistry === nothing && throw(ArgumentError("The model has no geochemical reactions"))
    return model.chemistry.selected[_time_index(model, t)]
end

@doc raw"""
    voltage(model, t) -> Float64

Potential difference ``\phi(0) - \phi(L)`` (V) between the left and the right boundary at
output time `t` (s).
"""
function voltage(model::Model1D, t::Real)
    φ = electric_potential(model, t)
    return φ[1] - φ[end]
end

@doc raw"""
    current_density(model, t) -> Vector{Float64}

Electric current density ``i = F \sum_i z_i J_i`` (A m⁻², per unit cross-sectional area of
the medium) in the ``+x`` direction, at output time `t` (s). The values refer to the midpoints
between neighbouring nodes, `(x[1:end-1] .+ x[2:end]) ./ 2`, and are computed from the
discrete fluxes of the transport equations.

With the [`ChargeConservation`](@ref) formulation, the current density is uniform. With the
[`Poisson`](@ref) formulation, the displacement current is not included; it is negligible
where the Debye length is not resolved by the grid.
"""
function current_density(model::Model1D, t::Real)
    u = model.results[_time_index(model, t)]
    d = model.problem.data
    iφ = d.nspecies + 1
    q = darcy_flux(d, u)
    i = zeros(Base.length(model.x) - 1)
    for k in eachindex(i)
        h = model.x[k + 1] - model.x[k]
        Δψ = d.f_RT * (u[iφ, k + 1] - u[iφ, k])
        for s in 1:(d.nspecies)
            i[k] += d.z[s] * sg_flux(d.Deff[s], d.z[s], u[s, k], u[s, k + 1], Δψ, q * h) / h
        end
    end
    return FARADAY .* i
end

@doc raw"""
    electroosmotic_flow(model, t) -> Float64

Electroosmotic water flux ``q`` (m s⁻¹, volume of water per unit cross-sectional area of the
medium and time) in the ``+x`` direction at output time `t` (s); see the section
Electroosmosis of [`Model1D`](@ref). The mean pore-water velocity is ``q / n``.
"""
electroosmotic_flow(model::Model1D, t::Real) = darcy_flux(model.problem.data, model.results[_time_index(model, t)])
