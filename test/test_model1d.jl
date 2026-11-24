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
end
