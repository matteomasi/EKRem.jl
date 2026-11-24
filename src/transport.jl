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
#     −∇·((ε_eff / F) ∇φ) = n Σᵢ zᵢ cᵢ          (Poisson equation divided by F)

"""
    PROGRESS_INTERVAL

Minimum wall-clock time (s) between two updates of the progress line shown during a run
when `progress=true`. The default is 0.5 s; 
It can be changed by setting the variable: 
`EKRem.PROGRESS_INTERVAL[] = 5.0`.
"""
const PROGRESS_INTERVAL = Ref(0.5)

"Parameters passed to the VoronoiFVM callbacks."
struct NPPData
    nspecies::Int
    porosity::Float64
    Deff::Vector{Float64}       # n τ D (m² s⁻¹)
    z::Vector{Float64}
    f_RT::Float64               # F / (R T) (V⁻¹)
    ε_F::Float64                # ε_eff / F (mol V⁻¹ m⁻¹)
    bvalues::Matrix{Float64}    # Dirichlet values (unknown × boundary region)
end

function npp_flux!(f, u, edge, d::NPPData)
    iφ = d.nspecies + 1
    Δψ = d.f_RT * (u[iφ, 2] - u[iφ, 1])
    for i in 1:(d.nspecies)
        f[i] = sg_flux(d.Deff[i], d.z[i], u[i, 1], u[i, 2], Δψ)
    end
    f[iφ] = d.ε_F * (u[iφ, 1] - u[iφ, 2])
    return nothing
end

function npp_storage!(f, u, node, d::NPPData)
    for i in 1:(d.nspecies)
        f[i] = d.porosity * u[i]
    end
    return nothing
end

function npp_reaction!(f, u, node, d::NPPData)
    f[d.nspecies + 1] = -d.porosity * sum(d.z[i] * u[i] for i in 1:(d.nspecies))
    return nothing
end

function npp_bcondition!(f, u, bnode, d::NPPData)
    ireg = bnode.region
    for i in 1:(d.nspecies + 1)
        boundary_dirichlet!(f, u, bnode, i, ireg, d.bvalues[i, ireg])
    end
    return nothing
end

"""
Nernst–Planck–Poisson problem on a grid: VoronoiFVM system, parameters and characteristic
scales of the unknowns (used for the error tolerances of the time integration).
"""
struct TransportProblem
    system::VoronoiFVM.System
    data::NPPData
    scales::Vector{Float64}
end

"""
    TransportProblem(grid, species, medium, ::Poisson, bvalues, scales; temperature)

Assemble the Nernst–Planck–Poisson system. `bvalues` is a (N+1) × (number of boundary
regions) matrix of Dirichlet values for the concentrations and the potential, and `scales`
contains the characteristic scale of each unknown.
"""
function TransportProblem(
    grid::ExtendableGrid,
    species::AbstractVector{Species},
    medium::PorousMedium,
    ::Poisson,
    bvalues::AbstractMatrix,
    scales::AbstractVector;
    temperature::Real,
)
    nspecies = length(species)
    data = NPPData(
        nspecies,
        medium.porosity,
        [effective_diffusion(medium, s) for s in species],
        [Float64(s.z) for s in species],
        FARADAY / (GAS_CONSTANT * temperature),
        effective_permittivity(medium) / FARADAY,
        Matrix{Float64}(bvalues),
    )
    system = VoronoiFVM.System(
        grid;
        species=collect(1:(nspecies + 1)),
        flux=npp_flux!,
        storage=npp_storage!,
        reaction=npp_reaction!,
        bcondition=npp_bcondition!,
        data=data,
        unknown_storage=:dense,
    )
    return TransportProblem(system, data, Vector{Float64}(scales))
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
    solve_transport(problem, u0, times; solver=:Rodas5P, reltol=1e-4, abstol=1e-6, progress=false)
        -> Vector{Matrix{Float64}}

Integrate the transport problem from the state `u0` (unknowns × nodes) at `times[1]` and
return the states at all `times` (the first one being `u0`).

The semi-discrete system is a differential-algebraic system in mass-matrix form (the
Poisson equation has no storage term). It is integrated by the method `solver` of
OrdinaryDiffEq.jl (one of [`SOLVERS`](@ref)), with adaptive time steps and error control;
the time steps end exactly at `times`. The error tolerances are `reltol` (relative) and
`abstol` times the largest concentration of each species (absolute).

The potential is excluded from the error control. It is an algebraic variable, determined
at each instant by the concentrations through the Poisson equation, so its accuracy
follows from theirs. Where the Debye length is not resolved, the potential responds to
charge on the time scale ε/σ (of the order of 1e-9 s), and including it in the error
estimate would force time steps of that size after any perturbation.

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
    state = VoronoiFVM.SystemState(problem.system)
    # The sparsity pattern of the Jacobian is detected at a state with positive
    # concentrations: where all concentrations vanish, the derivatives of the fluxes with
    # respect to the potential vanish too, and the pattern would be incomplete.
    jacval = Matrix{Float64}(u0)
    jacval[1:(nunknowns - 1), :] .= problem.scales[1:(nunknowns - 1)]
    f = SciMLBase.ODEFunction(state; jacval=vec(jacval), tjac=t[1])
    prob = SciMLBase.ODEProblem(f, vec(Matrix{Float64}(u0)), (t[1], t[end]), state)
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
    )
    if !SciMLBase.successful_retcode(sol)
        progress && report.last > 0 && println(report.io)
        error("The time integration stopped at t = $(sol.t[end]) s (return code $(sol.retcode))")
    end
    progress && show_progress!(report, sol.t[end]; final=true)
    return [reshape(u, nunknowns, nnodes) for u in sol.u]
end
