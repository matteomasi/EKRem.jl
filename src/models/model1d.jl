# One-dimensional model: columns and laboratory cells.

"""
    Model1D(; length, grid, species, initial, left, right, potential=Poisson(),
            porosity=1.0, tortuosity=1.0, relative_permittivity=WATER_RELATIVE_PERMITTIVITY,
            temperature=298.15)

One-dimensional model of electrokinetic transport in a water-saturated porous medium
between two boundaries at x = 0 (`left`) and x = `length` (`right`), in m.

# Arguments
- `grid`: node coordinates (a vector from 0 to `length`), or a named tuple
  `(h_min=…, h_max=…, growth=…)` for a grid refined toward both ends
  (see [`grid_refined`](@ref)).
- `species`: vector of [`Species`](@ref).
- `initial`: uniform initial pore-water concentrations (mol m⁻³), as a dictionary that maps
  species names to concentrations.
- `left`, `right`: boundary conditions, e.g. [`Reservoir`](@ref). The potential must be
  fixed at one boundary at least; a current density can be applied at the other
  one with the [`ChargeConservation`](@ref) formulation.
- `potential`: formulation of the potential equation: [`Poisson`](@ref),
  [`ChargeConservation`](@ref) or [`FixedField`](@ref).
- `porosity`, `tortuosity`, `relative_permittivity`: properties of the medium
  (see [`PorousMedium`](@ref)).
- `temperature`: temperature (K).

The initial and boundary compositions must be electroneutral. The model is run with
[`run!`](@ref), and the results are obtained with [`concentration`](@ref),
[`electric_potential`](@ref), [`voltage`](@ref) and [`current_density`](@ref) at the times
in [`output_times`](@ref).
"""
mutable struct Model1D
    x::Vector{Float64}
    species::Vector{Species}
    medium::PorousMedium
    temperature::Float64
    potential::PotentialModel
    initial::Vector{Float64}
    left::Reservoir
    right::Reservoir
    problem::TransportProblem
    state::Matrix{Float64}
    time::Float64
    times::Vector{Float64}
    results::Vector{Matrix{Float64}}
end

function Model1D(;
    length::Real,
    grid,
    species::AbstractVector{Species},
    initial::AbstractDict,
    left::Reservoir,
    right::Reservoir,
    potential::PotentialModel=Poisson(),
    porosity::Real=1.0,
    tortuosity::Real=1.0,
    relative_permittivity::Real=WATER_RELATIVE_PERMITTIVITY,
    temperature::Real=298.15,
)
    x = _coordinates(length, grid)
    names = [s.name for s in species]
    allunique(names) || throw(ArgumentError("Species names must be unique"))
    medium = PorousMedium(; porosity, tortuosity, relative_permittivity)
    _check_electric_conditions(potential, left, right)

    c0 = concentration_vector(species, initial)
    cl = concentration_vector(species, left.concentrations)
    cr = concentration_vector(species, right.concentrations)
    for (label, c) in (("initial", c0), ("left boundary", cl), ("right boundary", cr))
        imbalance = charge_imbalance(species, c)
        imbalance <= 1.0e-8 ||
            throw(ArgumentError("The $label composition is not electroneutral (relative charge imbalance $imbalance)"))
    end

    # Characteristic scales: largest concentration of each species (with a floor for species
    # absent from all compositions) and R T / F for the potential
    cmax = max.(c0, cl, cr)
    floor_c = 1.0e-9 * max(maximum(cmax), 1.0)
    scales = vcat(max.(cmax, floor_c), GAS_CONSTANT * temperature / FARADAY)

    electric(b) = fixed_potential(b) ? b.potential : b.current_density
    bvalues = hcat(vcat(cl, electric(left)), vcat(cr, electric(right)))
    problem = TransportProblem(
        simplexgrid(x),
        species,
        medium,
        potential,
        bvalues,
        [fixed_potential(left), fixed_potential(right)],
        scales;
        temperature=temperature,
    )

    # Initial state: uniform concentrations and boundary values at the boundary nodes. The
    # potential is computed from the potential equation, starting from a linear profile
    # (with an applied current density, the field is estimated from Ohm's law with the
    # initial conductivity).
    nspecies = Base.length(species)
    state = zeros(nspecies + 1, Base.length(x))
    state[1:nspecies, :] .= c0
    state[1:nspecies, 1] .= cl
    state[1:nspecies, end] .= cr
    state[end, :] .= _initial_potential(x, left, right, problem.data.Deff, [s.z for s in species], c0, temperature)
    state = solve_potential(problem, state)

    return Model1D(
        x, collect(species), medium, temperature, potential, c0, left, right, problem, state, 0.0, [0.0], [copy(state)]
    )
end

function _check_electric_conditions(potential::PotentialModel, left::Reservoir, right::Reservoir)
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

function _initial_potential(x, left::Reservoir, right::Reservoir, Deff, z, c0, temperature)
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
    compositions = (
        model.initial,
        concentration_vector(model.species, model.left.concentrations),
        concentration_vector(model.species, model.right.concentrations),
    )
    λ = maximum(
        debye_length(c, z; εr=model.medium.relative_permittivity, T=model.temperature) for
        c in compositions if any(>(0), abs.(z) .* c)
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
    electric(b) = fixed_potential(b) ? "φ = $(b.potential) V" : "i = $(b.current_density) A/m²"
    println(io, "  boundaries: left ", electric(m.left), ", right ", electric(m.right))
    return print(
        io, "  potential: ", nameof(typeof(m.potential)), "; λ_D,max / h_min = ", round(debye_ratio(m); sigdigits=3)
    )
end

"""
    run!(model, times; solver=:Rodas5P, reltol=1e-4, abstol=1e-6, progress=true) -> model

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
  that move by migration; one Jacobian factorisation per step.
- `:FBDF`: backward differentiation formulas of variable order (1 to 5). Fewer Jacobian
  evaluations per step; suited to long simulations with slowly varying solutions.

Solvers can be mixed for different parts of the simulation by calling `run!` multiple 
times with different `solver` arguments, e.g., see the example below.

# Example
```julia
run!(model, [3600.0, 86400.0])                  # Rodas5P, output after 1 h and 1 d
run!(model, [7 * 86400.0]; solver=:FBDF)        # continue for 7 d with FBDF
```
"""
function run!(
    model::Model1D, times; solver::Symbol=:Rodas5P, reltol::Real=1.0e-4, abstol::Real=1.0e-6, progress::Bool=true
)
    t = sort(unique(collect(Float64, times)))
    t[1] > model.time || throw(ArgumentError("Output times must be later than the current time $(model.time) s"))
    states = solve_transport(model.problem, model.state, vcat(model.time, t); solver, reltol, abstol, progress)
    for (k, tk) in enumerate(t)
        push!(model.times, tk)
        push!(model.results, copy(states[k + 1]))
    end
    model.state = copy(model.results[end])
    model.time = model.times[end]
    return model
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
    electric_potential(model, t) -> Vector{Float64}

Electric potential (V) at the nodes, at output time `t` (s).
"""
electric_potential(model::Model1D, t::Real) = model.results[_time_index(model, t)][end, :]

"""
    voltage(model, t) -> Float64

Potential difference φ(0) − φ(L) (V) between the left and the right boundary at output
time `t` (s).
"""
function voltage(model::Model1D, t::Real)
    φ = electric_potential(model, t)
    return φ[1] - φ[end]
end

"""
    current_density(model, t) -> Vector{Float64}

Electric current density i = F Σᵢ zᵢ Jᵢ (A m⁻², per unit cross-sectional area of the
medium) in the +x direction, at output time `t` (s). The values refer to the midpoints
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
    i = zeros(Base.length(model.x) - 1)
    for k in eachindex(i)
        h = model.x[k + 1] - model.x[k]
        Δψ = d.f_RT * (u[iφ, k + 1] - u[iφ, k])
        for s in 1:(d.nspecies)
            i[k] += d.z[s] * sg_flux(d.Deff[s], d.z[s], u[s, k], u[s, k + 1], Δψ) / h
        end
    end
    return FARADAY .* i
end
