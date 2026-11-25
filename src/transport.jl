# Nernst–Planck–Poisson system is assembled with VoronoiFVM.jl and the solution is computed
# over a time interval.
# The VoronoiFVM callbacks are defined on single edges and nodes and therefore do not depend on the
# spatial dimension.
#
# The unknowns, calculated at each node are:
# - the pore-water concentrations cᵢ (mol m⁻³) of the N species
# - electric potential φ (V), stored as unknowns 1...N and N+1.
#
# Equations, per unit bulk volume:
#     n ∂cᵢ/∂t + ∇·Jᵢ = 0,    Jᵢ = −Dᵢ* (∇cᵢ + zᵢ cᵢ ∇ψ),    Dᵢ* = n τ Dᵢ,  ψ = F φ / (R T)
# At a boundary with an electrode compartment (well-mixed solution of volume V per
# cross-sectional area A), the boundary node represents the compartment:
#     (V/A) ∂cᵢ/∂t + Jᵢ·ν + (Q/A) cᵢ = Rᵢ / A
# with ν the outward normal of the compartment toward the column, Q the flushing rate and Rᵢ
# the production by electrolysis.
# One of the following equations is solved for the potential:
#     −∇·((ε_eff / F) ∇φ) = n Σᵢ zᵢ cᵢ          Poisson equation divided by F
#     ∇·(Σᵢ zᵢ Jᵢ) = 0                         charge conservation divided by F
#     −∇·(∇φ) = 0                              fixed field (Laplace equation)
# The formulation is selected at run time by a field of the callback data, so that all
# formulations share the compiled code of the time integration.

"""
    PROGRESS_INTERVAL

Minimum wall-clock time (s) between two updates of the progress line shown during a run
when `progress=true`. The default is 0.5 s; 
It can be changed by setting the variable: 
`EKRem.PROGRESS_INTERVAL[] = 5.0`.
"""
const PROGRESS_INTERVAL = Ref(0.5)

# Codes of the formulations of the potential equation in NPPData
const POISSON = 1
const CHARGE_CONSERVATION = 2
const FIXED_FIELD = 3

formulation_code(::Poisson) = POISSON
formulation_code(::ChargeConservation) = CHARGE_CONSERVATION
formulation_code(::FixedField) = FIXED_FIELD

"Parameters passed to the VoronoiFVM callbacks."
struct NPPData
    nspecies::Int
    formulation::Int            # POISSON, CHARGE_CONSERVATION or FIXED_FIELD
    porosity::Float64
    Deff::Vector{Float64}       # n τ D (m² s⁻¹)
    z::Vector{Float64}
    f_RT::Float64               # F / (R T) (V⁻¹)
    ε_F::Float64                # ε_eff / F (mol V⁻¹ m⁻¹)
    bvalues::Matrix{Float64}    # boundary values (unknown × boundary region), see TransportProblem
    fixed_potential::Vector{Bool}   # per boundary region: potential (true) or current density (false)
    compartment::Vector{Bool}       # per boundary region: electrode compartment (true) or reservoir
    depth::Vector{Float64}          # compartment volume per cross-sectional area V/A (m)
    flushing::Vector{Float64}       # flushing rate per cross-sectional area Q/A (m s⁻¹)
    production::Matrix{Float64}     # production by electrolysis R/A (mol m⁻² s⁻¹), species × region
end

function npp_flux!(f, u, edge, d::NPPData)
    iφ = d.nspecies + 1
    Δψ = d.f_RT * (u[iφ, 2] - u[iφ, 1])
    charge_flux = zero(Δψ)
    for i in 1:(d.nspecies)
        f[i] = sg_flux(d.Deff[i], d.z[i], u[i, 1], u[i, 2], Δψ)
        charge_flux += d.z[i] * f[i]
    end
    if d.formulation == POISSON
        f[iφ] = d.ε_F * (u[iφ, 1] - u[iφ, 2])
    elseif d.formulation == CHARGE_CONSERVATION
        f[iφ] = charge_flux
    else
        f[iφ] = u[iφ, 1] - u[iφ, 2]
    end
    return nothing
end

function npp_storage!(f, u, node, d::NPPData)
    for i in 1:(d.nspecies)
        f[i] = d.porosity * u[i]
    end
    return nothing
end

function npp_reaction!(f, u, node, d::NPPData)
    if d.formulation == POISSON
        f[d.nspecies + 1] = -d.porosity * sum(d.z[i] * u[i] for i in 1:(d.nspecies))
    end
    return nothing
end

function npp_bstorage!(f, u, bnode, d::NPPData)
    ireg = bnode.region
    if d.compartment[ireg]
        for i in 1:(d.nspecies)
            f[i] = d.depth[ireg] * u[i]
        end
    end
    return nothing
end

npp_bcondition!(f, u, bnode, d::NPPData) = boundary_terms!(f, u, bnode.region, d)

# Terms of the boundary region `ireg`. Fixed values (Dirichlet conditions) are imposed by a
# penalty, as in VoronoiFVM.
function boundary_terms!(f, u, ireg, d::NPPData)
    for i in 1:(d.nspecies)
        if d.compartment[ireg]
            f[i] += d.flushing[ireg] * u[i] - d.production[i, ireg]
        else
            f[i] += PENALTY * (u[i] - d.bvalues[i, ireg])
        end
    end
    iφ = d.nspecies + 1
    if d.fixed_potential[ireg]
        f[iφ] += PENALTY * (u[iφ] - d.bvalues[iφ, ireg])
    else
        # Applied current density i: inflow i / F of the charge equation
        f[iφ] -= d.bvalues[iφ, ireg] / FARADAY
    end
    return nothing
end

const PENALTY = VoronoiFVM.Dirichlet(Float64)

"""
Nernst–Planck–Poisson problem on a grid: VoronoiFVM system, parameters, characteristic
scales of the unknowns (used for the error tolerances of the time integration), and node
coordinates of the one-dimensional grid.
"""
struct TransportProblem
    system::VoronoiFVM.System
    data::NPPData
    scales::Vector{Float64}
    x::Vector{Float64}
end

@doc raw"""
    TransportProblem(grid, species, medium, potential, bvalues, fixed_potential, scales;
                     temperature, depth, flushing, production)

Assemble the Nernst–Planck system with the equation `potential` for the electric
potential. `bvalues` is a ``(N+1) \times`` (number of boundary regions) matrix of boundary
values: Dirichlet values of the concentrations, and either the potential (V) or the current
density applied at the boundary (A m⁻²), as indicated by the vector `fixed_potential`.
`scales` contains the characteristic scale of each unknown.

Boundary regions with a positive `depth` (``V/A``, m) are electrode compartments, with the
flushing rate `flushing` (``Q/A``, m s⁻¹) and the production `production` (``R/A``,
mol m⁻² s⁻¹, species × region); their concentrations are not fixed.
"""
function TransportProblem(
    grid::ExtendableGrid,
    species::AbstractVector{Species},
    medium::PorousMedium,
    potential::PotentialModel,
    bvalues::AbstractMatrix,
    fixed_potential::AbstractVector{Bool},
    scales::AbstractVector;
    temperature::Real,
    depth::AbstractVector=zeros(size(bvalues, 2)),
    flushing::AbstractVector=zeros(size(bvalues, 2)),
    production::AbstractMatrix=zeros(length(species), size(bvalues, 2)),
)
    nspecies = length(species)
    data = NPPData(
        nspecies,
        formulation_code(potential),
        medium.porosity,
        [effective_diffusion(medium, s) for s in species],
        [Float64(s.z) for s in species],
        FARADAY / (GAS_CONSTANT * temperature),
        effective_permittivity(medium) / FARADAY,
        Matrix{Float64}(bvalues),
        Vector{Bool}(fixed_potential),
        depth .> 0,
        Vector{Float64}(depth),
        Vector{Float64}(flushing),
        Matrix{Float64}(production),
    )
    system = VoronoiFVM.System(
        grid;
        species=collect(1:(nspecies + 1)),
        flux=npp_flux!,
        storage=npp_storage!,
        reaction=npp_reaction!,
        bcondition=npp_bcondition!,
        bstorage=npp_bstorage!,
        data=data,
        unknown_storage=:dense,
    )
    x = vec(grid[ExtendableGrids.Coordinates])
    return TransportProblem(system, data, Vector{Float64}(scales), x)
end

# Progress report
mutable struct Progress
    t0::Float64
    t1::Float64
    start::Float64
    last::Float64
    interval::Float64
    io::IO
end

Progress(t0, t1; io::IO=stdout, interval::Real=PROGRESS_INTERVAL[]) = Progress(t0, t1, time(), 0.0, interval, io)

function show_progress!(p::Progress, t; final=false)
    elapsed = time() - p.start
    final || elapsed - p.last >= p.interval || return nothing
    percent = 100 * (t - p.t0) / (p.t1 - p.t0)
    message = @sprintf("  Simulated %3.0f %% (t = %.4g s of %.4g s), %.1f s elapsed", percent, t, p.t1, elapsed)
    print(p.io, "\r", rpad(message, 72))
    final && println(p.io)
    flush(p.io)
    p.last = elapsed
    return nothing
end

function progress_callback(p::Progress)
    function affect!(integrator)
        show_progress!(p, integrator.t)
        SciMLBase.u_modified!(integrator, false)
        return nothing
    end
    return SciMLBase.DiscreteCallback((u, t, integrator) -> true, affect!; save_positions=(false, false))
end

"Time integration methods available for the transport equations."
const SOLVERS = (:Rodas5P, :FBDF)

# UMFPACK is selected explicitly: the default sparse factorisation of LinearSolve.jl requires
# a fixed sparsity pattern, which the VoronoiFVM matrix does not guarantee.
function time_integrator(solver::Symbol)
    solver === :Rodas5P && return Rodas5P(; linsolve=UMFPACKFactorization())
    solver === :FBDF && return FBDF(; linsolve=UMFPACKFactorization())
    return throw(ArgumentError("Unknown solver :$solver; the available solvers are $(join(repr.(SOLVERS), ", "))"))
end

"""
    ode_function(problem, u0, t0) -> (state, f)

VoronoiFVM system state and ODE function (right-hand side, Jacobian and mass matrix) of
the semi-discrete transport problem.
"""
function ode_function(problem::TransportProblem, u0::AbstractMatrix, t0)
    nunknowns = size(u0, 1)
    state = VoronoiFVM.SystemState(problem.system)
    # The sparsity pattern of the Jacobian is detected at a state with positive
    # concentrations: where all concentrations vanish, the derivatives of the fluxes with
    # respect to the potential vanish too, and the pattern would be incomplete.
    jacval = Matrix{Float64}(u0)
    jacval[1:(nunknowns - 1), :] .= problem.scales[1:(nunknowns - 1)]
    vfvm = SciMLBase.ODEFunction(state; jacval=vec(jacval), tjac=t0)
    # VoronoiFVM evaluates the Jacobian with every right-hand side; the Rosenbrock and BDF
    # methods need the right-hand side more often than the Jacobian, so the residual is
    # evaluated separately (see `residual!`).
    workspace = ResidualWorkspace(nunknowns)
    rhs!(du, u, p, t) = (residual!(du, u, problem, workspace); du .*= -1; nothing)
    f = SciMLBase.ODEFunction(rhs!; jac=vfvm.jac, jac_prototype=vfvm.jac_prototype, mass_matrix=vfvm.mass_matrix)
    return state, f
end

struct ResidualWorkspace
    u::Matrix{Float64}      # unknowns at the two nodes of an edge
    f::Vector{Float64}
end

ResidualWorkspace(n) = ResidualWorkspace(zeros(n, 2), zeros(n))

@doc raw"""
    residual!(r, u, problem, workspace)

Residual of the semi-discrete equations without the storage terms, as assembled by
VoronoiFVM, for the state vector `u` (unknowns × nodes, column-major), on the
one-dimensional grid of the problem: edge fluxes times ``1/h``, reactions times the Voronoi
volumes of the nodes, and boundary terms. Unlike VoronoiFVM, no derivatives are computed.
"""
function residual!(r, u, problem::TransportProblem, w::ResidualWorkspace)
    d, x = problem.data, problem.x
    n = size(w.u, 1)
    U = reshape(u, n, :)
    R = reshape(r, n, :)
    fill!(R, 0.0)
    for k in 1:(length(x) - 1)
        h = x[k + 1] - x[k]
        @views w.u[:, 1] .= U[:, k]
        @views w.u[:, 2] .= U[:, k + 1]
        fill!(w.f, 0.0)
        npp_flux!(w.f, w.u, nothing, d)
        @views R[:, k] .+= w.f ./ h
        @views R[:, k + 1] .-= w.f ./ h
        # Reaction terms over the two halves of the edge (Voronoi volumes)
        for j in (k, k + 1)
            fill!(w.f, 0.0)
            @views npp_reaction!(w.f, U[:, j], nothing, d)
            @views R[:, j] .+= w.f .* (h / 2)
        end
    end
    for (j, ireg) in ((1, 1), (length(x), 2))
        fill!(w.f, 0.0)
        @views boundary_terms!(w.f, U[:, j], ireg, d)
        @views R[:, j] .+= w.f
    end
    return r
end

"""
    solve_potential(problem, u0, t0=0.0) -> Matrix{Float64}

Copy of the state `u0` (unknowns × nodes) with the potential computed by
[`solve_potential!`](@ref).
"""
function solve_potential(problem::TransportProblem, u0::AbstractMatrix, t0=0.0)
    state, f = ode_function(problem, u0, t0)
    u = solve_potential!(vec(Matrix{Float64}(u0)), problem, f, state, t0)
    return reshape(u, size(u0))
end

"""
    solve_potential!(u, problem, f, state, t) -> u

Solve the potential equation for the concentrations in the state vector `u` (Newton
iteration), and replace the potential in `u` by the solution. Not applied with the
Poisson formulation, where the potential is left as computed by the time integration.
"""
function solve_potential!(u, problem::TransportProblem, f, state, t; maxiters=50)
    # With the Poisson formulation, the potential is set by charge imbalances far below the
    # tolerances, which this solution would amplify.
    problem.data.formulation == POISSON && return u
    nunknowns = problem.data.nspecies + 1
    scale = problem.scales[end]
    iφ = nunknowns:nunknowns:length(u)
    rhs = similar(u)
    jac = copy(f.jac_prototype)
    for _ in 1:maxiters
        f.f(rhs, u, state, t)       # rhs = −residual, which vanishes in the rows of the potential
        f.jac(jac, u, state, t)
        δ = jac[iφ, iφ] \ rhs[iφ]
        # Damping limits the change of the potential to 10 thermal voltages per iteration
        δmax = maximum(abs, δ; init=0.0)
        damping = min(1.0, 10 * scale / δmax)
        u[iφ] .-= damping .* δ
        damping == 1 && δmax <= 1.0e-10 * scale && return u
    end
    return error("The potential equation could not be solved at t = $t s (no convergence in $maxiters iterations)")
end

@doc raw"""
    solve_transport(problem, u0, times; solver=:Rodas5P, reltol=1e-4, abstol=1e-6, progress=false)
        -> Vector{Matrix{Float64}}

Integrate the transport problem from the state `u0` (unknowns × nodes) at `times[1]` and
return the states at all `times` (the first one being `u0`).

The semi-discrete system is a differential-algebraic system in mass-matrix form (the
equation of the potential has no storage term). It is integrated by the method `solver` of
OrdinaryDiffEq.jl (one of [`SOLVERS`](@ref)), with adaptive time steps and error control;
the time steps end exactly at `times`. The error tolerances are `reltol` (relative) and
`abstol` times the largest concentration of each species (absolute).

The potential is excluded from the error control. It is an algebraic variable, determined
at each instant by the concentrations, so its accuracy follows from theirs. With the
Poisson formulation, where the Debye length is not resolved, the potential responds to
charge on the time scale ``\varepsilon / \sigma`` (of the order of ``10^{-9}`` s), and including it in the error
estimate would force time steps of that size after any perturbation. With the other
formulations, the time integration solves the potential equation only approximately, and
the current density is sensitive to the remaining error; the potential of the initial and
the returned states is therefore computed by [`solve_potential!`](@ref).

The method does not guarantee non-negative concentrations. Undershoots smaller than the
absolute tolerance are accepted; time steps that produce larger negative concentrations
are rejected and repeated with a smaller step.

With `progress=true`, the progress and elapsed time are printed at intervals of
[`PROGRESS_INTERVAL`](@ref) seconds.
"""
function solve_transport(
    problem::TransportProblem,
    u0::AbstractMatrix,
    times;
    solver::Symbol=:Rodas5P,
    reltol=1.0e-4,
    abstol=1.0e-6,
    progress::Bool=false,
)
    algorithm = time_integrator(solver)
    t = collect(Float64, times)
    nunknowns, nnodes = size(u0)
    state, f = ode_function(problem, u0, t[1])
    u_start = solve_potential!(vec(Matrix{Float64}(u0)), problem, f, state, t[1])
    prob = SciMLBase.ODEProblem(f, u_start, (t[1], t[end]), state)
    atol = repeat(problem.scales .* abstol, nnodes)
    atol[nunknowns:nunknowns:end] .= Inf   # potential excluded from the error control
    # Lower bound of each unknown: −abstol for the concentrations, none for the potential
    lower = [mod1(k, nunknowns) < nunknowns ? -atol[k] : -Inf for k in eachindex(atol)]
    outofdomain(u, p, t) = any(u[k] < lower[k] for k in eachindex(u))
    report = Progress(t[1], t[end])
    sol = SciMLBase.solve(
        prob,
        algorithm;
        reltol=reltol,
        abstol=atol,
        isoutofdomain=outofdomain,
        saveat=t,
        tstops=t,
        callback=progress ? progress_callback(report) : nothing,
        # The initial potential is computed above; the default check would compare
        # the residuals of the Dirichlet rows (penalty method) with the tolerances.
        initializealg=SciMLBase.NoInit(),
    )
    if !SciMLBase.successful_retcode(sol)
        progress && report.last > 0 && println(report.io)
        error("The time integration stopped at t = $(sol.t[end]) s (return code $(sol.retcode))")
    end
    progress && show_progress!(report, sol.t[end]; final=true)
    return [reshape(solve_potential!(u, problem, f, state, tk), nunknowns, nnodes) for (u, tk) in zip(sol.u, sol.t)]
end
