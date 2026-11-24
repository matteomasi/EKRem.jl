# Tests of the electrode compartments: mass balances in a closed system, with and without
# geochemical reactions.

using Test
using EKRem

const WATER_SPECIES = [
    Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), Species("H+", 1, 9.31e-9), Species("OH-", -1, 5.27e-9)
]

@testset "Electrode compartments" begin
    nacl = Dict("Na+" => 10.0, "Cl-" => 10.0)
    i = 10.0   # A/m², applied at the left boundary (anode)
    column(left, right; area=1.0e-3) = Model1D(;
        length=0.02,
        grid=collect(range(0, 0.02; length=21)),
        species=WATER_SPECIES,
        initial=nacl,
        left,
        right,
        potential=ChargeConservation(),
        porosity=0.4,
        tortuosity=0.5,
        area,
    )

    @testset "Mass balances" begin
        # Closed compartments: Na+ and Cl- are conserved, and electrolysis adds i t / F of
        # H+ (anode) and of OH- (cathode); the compartments remain electroneutral
        model = column(
            Compartment(nacl; volume=1.0e-4, current_density=i), Compartment(nacl; volume=1.0e-4, potential=0.0)
        )
        run!(model, [600.0, 1800.0]; progress=false)
        for t in (600.0, 1800.0)
            Δ(name) = amount(model, name, t) - amount(model, name, 0.0)
            @test abs(Δ("Na+")) < 1.0e-12 * amount(model, "Na+", 0.0)
            @test abs(Δ("Cl-")) < 1.0e-12 * amount(model, "Cl-", 0.0)
            @test Δ("H+") ≈ i * t / FARADAY rtol = 1.0e-10
            @test Δ("OH-") ≈ i * t / FARADAY rtol = 1.0e-10
            c = [concentration(model, s.name, t) for s in WATER_SPECIES]
            charge = c[1] .- c[2] .+ c[3] .- c[4]
            @test abs(charge[1]) < 1.0e-10 && abs(charge[end]) < 1.0e-10
            @test all(current_density(model, t) .≈ i)
        end
        # H+ accumulates in the anode compartment, OH- in the cathode compartment
        @test concentration(model, "H+", 1800.0)[1] > 0.1 * i * 1800.0 / FARADAY / (1.0e-4 / 1.0e-3)
        @test concentration(model, "OH-", 1800.0)[end] > concentration(model, "OH-", 1800.0)[1]

        # Current applied at the cathode (negative current density at the right boundary):
        # OH- is produced at the right boundary
        cathode = column(Reservoir(nacl; potential=0.0), Compartment(nacl; volume=1.0e-4, current_density=-i))
        run!(cathode, [600.0]; progress=false)
        @test amount(cathode, "OH-", 600.0) - amount(cathode, "OH-", 0.0) > 0.5 * i * 600.0 / FARADAY

        # Flushing with deionised water removes Na+ from the system
        flushed = column(
            Compartment(nacl; volume=1.0e-4, flow=1.0e-7, current_density=i),
            Compartment(nacl; volume=1.0e-4, flow=1.0e-7, potential=0.0),
        )
        run!(flushed, [600.0]; progress=false)
        @test amount(flushed, "Na+", 600.0) < amount(flushed, "Na+", 0.0)
    end

    @testset "Input checks" begin
        compartment = Compartment(nacl; volume=1.0e-4, current_density=i)
        reservoir = Reservoir(nacl; potential=0.0)
        @test_throws ArgumentError column(compartment, reservoir; area=nothing)              # area required
        @test_throws ArgumentError column(                                                   # voltage control
            Compartment(nacl; volume=1.0e-4, potential=1.0),
            reservoir,
        )
        @test_throws ArgumentError Compartment(nacl; volume=0.0, current_density=i)
        @test_throws ArgumentError Compartment(nacl; volume=1.0e-4, flow=-1.0, current_density=i)
        @test_throws ArgumentError Compartment(nacl; volume=1.0e-4)                          # no electric condition
        @test_throws ArgumentError Model1D(;                                                 # no H+ for electrolysis
            length=0.02,
            grid=collect(range(0, 0.02; length=21)),
            species=WATER_SPECIES[[1, 2, 4]],
            initial=nacl,
            left=compartment,
            right=reservoir,
            potential=ChargeConservation(),
            area=1.0e-3,
        )
    end

    lib = phreeqcrm_library()
    if lib === nothing
        @test_skip false
    else
        @testset "Compartments with geochemical reactions" begin
            database = get(ENV, "EKREM_PHREEQC_DATABASE", joinpath(dirname(dirname(lib)), "database", "phreeqc.dat"))
            # Calcite in the column; the catholyte is kept at pH 3 by addition of HNO3
            input = """
            PHASES
            Fix_H+
                H+ = H+
                log_k 0
            SOLUTION 1
                pH 7 charge
                Na 10
                Cl 10
            EQUILIBRIUM_PHASES 1
                Calcite 0 0.1
            SOLUTION 2
                pH 7 charge
                Na 1
                Cl 1
            EQUILIBRIUM_PHASES 2
                Fix_H+ -3 HNO3 10
            SELECTED_OUTPUT 1
                -reset false
                -pH true
            END
            """
            model = Model1D(;
                length=0.02,
                grid=collect(range(0, 0.02; length=21)),
                initial=(solution=1, equilibrium_phases=1),
                left=Compartment(; initial=(solution=2,), volume=1.0e-4, current_density=i),
                right=Compartment(; initial=(solution=2, equilibrium_phases=2), volume=1.0e-4, potential=0.0),
                database,
                input,
                potential=ChargeConservation(),
                porosity=0.4,
                tortuosity=0.5,
                area=1.0e-3,
            )
            run!(model, [300.0]; coupling_step=60.0, progress=false)
            pH = selected_output(model, 300.0)["pH"]
            @test pH[1] < 4                                     # acidified anolyte
            @test abs(pH[end] - 3) < 1.0e-6                     # catholyte at pH 3
            @test concentration(model, "NO3-", 300.0)[end] > 0  # nitric acid added
            # Cl is conserved (Rodas5P conserves the amounts to round-off)
            cl(t) = amount(model, "Cl-", t) + amount(model, "HCl", t)
            @test cl(300.0) ≈ cl(0.0) rtol = 1.0e-9
            @test all(current_density(model, 300.0) .≈ i)
        end
    end
end
