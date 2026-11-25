# Unit tests of the one-dimensional model on a small grid. Complete simulations compared with
# analytical solutions are in the examples/ folder.

using Test
using EKRem

const SPECIES = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), Species("OH-", -1, 5.27e-9)]

# NaCl in the column; NaOH (no Cl-) in the left reservoir at 1 V; NaCl in the right reservoir
function small_model(; kwargs...)
    return Model1D(;
        length=0.01,
        grid=collect(range(0, 0.01; length=21)),
        species=SPECIES,
        initial=Dict("Na+" => 10.0, "Cl-" => 10.0),
        left=Reservoir(Dict("Na+" => 10.0, "OH-" => 10.0); potential=1.0),
        right=Reservoir(Dict("Na+" => 100.0, "Cl-" => 100.0); potential=0.0),
        kwargs...,
    )
end

@testset "Model1D" begin
    model = small_model()
    @test coordinates(model) == collect(range(0, 0.01; length=21))
    @test output_times(model) == [0.0]
    @test concentration(model, "Na+", 0.0)[2] == 10.0
    @test electric_potential(model, 0.0) ≈ range(1.0, 0.0; length=21)
    @test debye_ratio(model) < 1.0e-4
    @test_throws ArgumentError concentration(model, "K+", 0.0)
    @test_throws ArgumentError concentration(model, "Na+", 5.0)

    # Output times are sorted; repeated calls continue from the last state
    run!(model, [100.0, 10.0]; progress=false)
    run!(model, [1000.0]; progress=false)
    @test output_times(model) == [0.0, 10.0, 100.0, 1000.0]
    @test_throws ArgumentError run!(model, [500.0])
    @test_throws ArgumentError run!(model, [2000.0]; solver=:RK4)   # unknown solver

    na, cl, oh = (concentration(model, name, 1000.0) for name in ("Na+", "Cl-", "OH-"))
    # Electroneutrality: the Debye length (nm) is much smaller than the grid spacing (0.5 mm)
    @test maximum(abs.(na .- cl .- oh) ./ na) < 1.0e-8
    # Cl- is depleted at the left boundary: concentrations remain non-negative up to the
    # absolute tolerance of the time integration (1e-6 of the largest concentration)
    @test minimum(cl) >= -1.0e-6 * 100
    # The solution does not depend on the permittivity (εr = 1 versus water)
    vacuum = small_model(; relative_permittivity=1.0)
    run!(vacuum, [1000.0]; progress=false)
    @test concentration(vacuum, "Na+", 1000.0) ≈ na rtol = 1.0e-3
    # Both time integration methods give the same solution within their tolerance
    fbdf = small_model()
    run!(fbdf, [1000.0]; solver=:FBDF, progress=false)
    @test concentration(fbdf, "Na+", 1000.0) ≈ na rtol = 1.0e-3
end

@testset "Formulations of the potential equation" begin
    poisson = small_model()
    run!(poisson, [1000.0]; progress=false)

    # Charge conservation: electroneutral limit of the Poisson formulation
    cc = small_model(; potential=ChargeConservation())
    run!(cc, [1000.0]; progress=false)
    na, cl, oh = (concentration(cc, name, 1000.0) for name in ("Na+", "Cl-", "OH-"))
    @test maximum(abs.(na .- cl .- oh) ./ na) < 1.0e-12
    @test na ≈ concentration(poisson, "Na+", 1000.0) rtol = 1.0e-6
    @test electric_potential(cc, 1000.0) ≈ electric_potential(poisson, 1000.0) atol = 1.0e-6
    # The current density is uniform (∇·i = 0), at all output times
    for t in (0.0, 1000.0)
        i = current_density(cc, t)
        @test length(i) == length(coordinates(cc)) - 1
        @test all(i .≈ i[1])
    end
    @test current_density(cc, 1000.0)[1] ≈ current_density(poisson, 1000.0)[1] rtol = 1.0e-6

    # Fixed field: the potential stays linear, whatever the concentrations
    pf = small_model(; potential=FixedField())
    run!(pf, [1000.0]; progress=false)
    @test electric_potential(pf, 1000.0) ≈ range(1.0, 0.0; length=21) atol = 1.0e-12
    @test voltage(pf, 1000.0) ≈ 1.0
end

@testset "Ohm's law in a uniform electrolyte" begin
    # Uniform NaCl in the column and in both reservoirs: no concentration gradients, so
    # i = σ U / L with the bulk conductivity σ of the medium
    electrolyte = Dict("Na+" => 10.0, "Cl-" => 10.0)
    n, τ, L, U = 0.4, 0.5, 0.01, 0.5
    σ = conductivity([10.0, 10.0], [1, -1], n * τ .* [1.33e-9, 2.03e-9])
    uniform(potential, left) = Model1D(;
        length=L,
        grid=collect(range(0, L; length=11)),
        species=SPECIES,
        initial=electrolyte,
        left,
        right=Reservoir(electrolyte; potential=0.0),
        potential,
        porosity=n,
        tortuosity=τ,
    )
    for potential in (Poisson(), ChargeConservation())
        model = uniform(potential, Reservoir(electrolyte; potential=U))
        run!(model, [100.0]; progress=false)
        @test all(current_density(model, 100.0) .≈ σ * U / L)
    end
    # Current control: the voltage follows from the applied current density
    model = uniform(ChargeConservation(), Reservoir(electrolyte; current_density=σ * U / L))
    @test voltage(model, 0.0) ≈ U
    run!(model, [100.0]; progress=false)
    @test voltage(model, 100.0) ≈ U
end

@testset "Current control" begin
    # Current density of the voltage-controlled model at t = 0, applied at the left boundary
    i0 = current_density(small_model(; potential=ChargeConservation()), 0.0)[1]
    model = small_model(;
        potential=ChargeConservation(), left=Reservoir(Dict("Na+" => 10.0, "OH-" => 10.0); current_density=i0)
    )
    @test voltage(model, 0.0) ≈ 1.0
    run!(model, [1000.0]; progress=false)
    @test all(current_density(model, 1000.0) .≈ i0)
    @test electric_potential(model, 1000.0)[end] ≈ 0.0 atol = 1.0e-12
    @test voltage(model, 1000.0) > 0
end

@testset "Electroosmosis" begin
    # Trace species between two reservoirs (c₀ at x = 0, none at x = L) under a uniform field
    # (fixed-field formulation) and the electroosmotic water flux q = k_eo U / L. Each species
    # moves with the velocity v = q/n + τ D z F E / (R T), and its steady state is
    #     c(x) = c₀ (e^Pe − e^(Pe x/L)) / (e^Pe − 1),    Pe = v L / (τ D),
    # which the Scharfetter–Gummel flux reproduces at the nodes (constant velocity).
    L, U, n, τ, k_eo, c0 = 0.01, 0.2, 0.4, 0.5, 2.0e-8, 0.01
    trace = [Species("N", 0, 1.0e-9), Species("T+", 1, 1.0e-9), Species("A-", -1, 2.0e-9)]
    background = Dict("Na+" => 100.0, "Cl-" => 100.0)
    x = collect(range(0, L; length=21))
    model = Model1D(;
        length=L,
        grid=x,
        species=[SPECIES[1:2]; trace],
        initial=background,
        left=Reservoir(merge(background, Dict(s.name => c0 for s in trace)); potential=U),
        right=Reservoir(background; potential=0.0),
        potential=FixedField(),
        porosity=n,
        tortuosity=τ,
        electroosmotic_permeability=k_eo,
    )
    q = k_eo * U / L
    @test electroosmotic_flow(model, 0.0) ≈ q
    run!(model, [1.0e6]; progress=false)
    @test electroosmotic_flow(model, 1.0e6) ≈ q
    for s in trace
        v = q / n + τ * s.D * s.z * FARADAY / (GAS_CONSTANT * 298.15) * U / L
        Pe = v * L / (τ * s.D)
        exact = c0 .* (exp(Pe) .- exp.(Pe .* x ./ L)) ./ (exp(Pe) - 1)
        @test concentration(model, s.name, 1.0e6) ≈ exact atol = 1.0e-6 * c0
    end
    # Without electroosmotic permeability, there is no water flux
    @test electroosmotic_flow(small_model(), 0.0) == 0
    @test_throws ArgumentError small_model(; electroosmotic_permeability=Inf)
end

@testset "Progress report" begin
    io = IOBuffer()
    p = EKRem.Progress(0.0, 1000.0; io, interval=2.0)
    EKRem.show_progress!(p, 100.0)
    @test isempty(take!(io))                 # no update during the first interval
    p.start -= 3.0                           # as if the run had started 3 s ago
    EKRem.show_progress!(p, 500.0)
    @test occursin("50 % (t = 500 s of 1000 s)", String(take!(io)))
    EKRem.show_progress!(p, 1000.0; final=true)
    @test occursin("100 %", String(take!(io)))   # the final line is always printed

    # The default interval is taken from the global setting
    default = EKRem.PROGRESS_INTERVAL[]
    try
        EKRem.PROGRESS_INTERVAL[] = 10.0
        @test EKRem.Progress(0.0, 1.0).interval == 10.0
    finally
        EKRem.PROGRESS_INTERVAL[] = default
    end
end

@testset "Model1D input checks" begin
    grid = collect(range(0, 0.01; length=11))
    electrolyte = Dict("Na+" => 10.0, "Cl-" => 10.0)
    build(; length=0.01, initial=electrolyte) = Model1D(;
        length,
        grid,
        species=SPECIES,
        initial,
        left=Reservoir(electrolyte; potential=1.0),
        right=Reservoir(electrolyte; potential=0.0),
    )
    @test build() isa Model1D
    @test_throws ArgumentError build(; initial=Dict("Na+" => 10.0, "Cl-" => 9.0))    # not electroneutral
    @test_throws ArgumentError build(; initial=Dict("Na+" => 10.0, "Br-" => 10.0))   # unknown species
    @test_throws ArgumentError build(; length=0.02)                                   # grid does not span [0, L]

    # Electric boundary conditions
    @test_throws ArgumentError Reservoir(electrolyte)                                   # neither potential nor current
    @test_throws ArgumentError Reservoir(electrolyte; potential=1.0, current_density=1.0)
    current = Reservoir(electrolyte; current_density=1.0)
    controlled(potential; right=Reservoir(electrolyte; potential=0.0)) =
        Model1D(; length=0.01, grid, species=SPECIES, initial=electrolyte, left=current, right, potential)
    @test controlled(ChargeConservation()) isa Model1D
    @test_throws ArgumentError controlled(Poisson())                                    # requires charge conservation
    @test_throws ArgumentError controlled(FixedField())
    @test_throws ArgumentError controlled(ChargeConservation(); right=current)          # no fixed potential
end

@testset "Residual of the transport equations" begin
    # The residual evaluated by EKRem (without derivatives) agrees with the assembly of
    # VoronoiFVM, for the three formulations, with electroosmosis and for electrode
    # compartments
    c = Dict("Na+" => 10.0, "Cl-" => 10.0, "OH-" => 1.0e-3)
    species = [SPECIES; Species("H+", 1, 9.31e-9)]
    c["H+"] = 1.0e-3
    models = [
        small_model(),
        small_model(; potential=FixedField()),
        small_model(; electroosmotic_permeability=1.0e-8),
        Model1D(;
            length=0.01,
            grid=collect(range(0, 0.01; length=21)),
            species,
            initial=c,
            left=Compartment(c; volume=1.0e-4, flow=1.0e-7, current_density=5.0),
            right=Compartment(c; volume=1.0e-4, potential=0.0),
            potential=ChargeConservation(),
            area=1.0e-3,
            electroosmotic_permeability=-1.0e-8,   # water flux toward the left compartment
        ),
    ]
    for model in models
        u = copy(model.state)
        u[1:(end - 1), :] .*= 1 .+ 0.1 .* sin.((1:(size(u, 1) - 1)) .+ 7 .* (1:size(u, 2))')   # perturbed state
        u[end, :] .+= 0.01 .* cos.(1:size(u, 2))
        state, f = EKRem.ode_function(model.problem, u, 0.0)
        r, r_vfvm = similar(vec(u)), similar(vec(u))
        f.f(r, vec(u), state, 0.0)
        EKRem.VoronoiFVM.eval_rhs!(r_vfvm, vec(u), state, 0.0)
        @test r ≈ r_vfvm rtol = 1.0e-12
    end
end
