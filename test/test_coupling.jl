# Tests of the coupling of transport and geochemical reactions in Model1D, on a small grid.

using Test
using EKRem

const NACL_SYSTEM = """
SOLUTION 1
    pH 7 charge
    Na 10
    Cl 10
SOLUTION 2
    pH 7 charge
    Na 100
    Cl 100
END
"""

# Column with calcite, acid at the anode (left), NaCl at the cathode
const ACID_CALCITE_SYSTEM = """
SOLUTION 1
    pH 7 charge
    Na 10
    Cl 10
EQUILIBRIUM_PHASES 1
    Calcite 0 0.05
SOLUTION 2
    pH 2
    Na 10
    Cl 10 charge
SOLUTION 3
    pH 7 charge
    Na 10
    Cl 10
SELECTED_OUTPUT 1
    -reset false
    -pH true
    -equilibrium_phases Calcite
    -saturation_indices Calcite
END
"""

@testset "Coupling of transport and reactions" begin
    lib = phreeqcrm_library()
    if lib === nothing
        @info "No PhreeqcRM library set (see `set_phreeqcrm_library`): coupling tests skipped"
        @test_skip false
    else
        database = get(ENV, "EKREM_PHREEQC_DATABASE", joinpath(dirname(dirname(lib)), "database", "phreeqc.dat"))
        grid = collect(range(0, 0.01; length=21))
        column(input, initial, left, right; kwargs...) = Model1D(;
            length=0.01,
            grid,
            initial,
            left=Reservoir(; solution=left, potential=0.5),
            right=Reservoir(; solution=right, potential=0.0),
            database,
            input,
            potential=ChargeConservation(),
            porosity=0.4,
            tortuosity=0.5,
            kwargs...,
        )

        @testset "Non-reactive solution" begin
            # NaCl only: the coupled model agrees with transport alone (Na+ and Cl-, with the
            # concentrations and diffusion coefficients of the coupled model)
            coupled = column(NACL_SYSTEM, (solution=1,), 2, 1)
            @test "H2O" ∉ [s.name for s in coupled.species]
            @test isempty(selected_output(coupled, 0.0))                       # no SELECTED_OUTPUT block
            run!(coupled, [1000.0]; coupling_step=100.0, progress=false)
            species = filter(s -> s.name in ("Na+", "Cl-"), coupled.species)
            nacl(k) = (c=concentration(coupled, "Na+", 0.0)[k]; Dict("Na+" => c, "Cl-" => c))
            transport = Model1D(;
                length=0.01,
                grid,
                species,
                initial=nacl(2),
                left=Reservoir(nacl(1); potential=0.5),
                right=Reservoir(nacl(21); potential=0.0),
                potential=ChargeConservation(),
                porosity=0.4,
                tortuosity=0.5,
            )
            run!(transport, [1000.0]; progress=false)
            for name in ("Na+", "Cl-")
                @test concentration(coupled, name, 1000.0) ≈ concentration(transport, name, 1000.0) rtol = 1.0e-4
            end

            # A subset of species, with user-defined diffusion coefficients
            subset = column(NACL_SYSTEM, (solution=1,), 2, 1; species=[Species("Na+", 1, 1.0e-9), "Cl-", "H+", "OH-"])
            @test [s.name for s in subset.species] == ["Na+", "Cl-", "H+", "OH-"]
            @test subset.species[1].D == 1.0e-9
        end

        @testset "Dissolution of calcite by acid" begin
            model = column(ACID_CALCITE_SYSTEM, (solution=1, equilibrium_phases=1), 2, 3)
            run!(model, [1000.0, 2000.0]; coupling_step=100.0, solver=:FBDF, progress=false)
            so = selected_output(model, 2000.0)
            calcite, si, pH = so["Calcite"], so["si_Calcite"], so["pH"]
            interior = 2:20
            # Calcite dissolves near the anode, where the acid enters, and remains elsewhere
            # in equilibrium with the pore water
            @test calcite[2] < 0.05
            @test all(abs.(si[interior][calcite[interior] .> 0]) .< 1.0e-6)
            @test pH[2] < pH[11]
            # Electroneutrality, non-negative concentrations and uniform current
            z = [s.z for s in model.species]
            c = model.results[end][1:length(z), :]
            @test all(abs.(z' * c) .<= 1.0e-6 .* (abs.(z)' * c))
            @test all(c .>= 0)
            i = current_density(model, 2000.0)
            @test all(i .≈ i[1])
        end

        @testset "Input checks" begin
            model = column(NACL_SYSTEM, (solution=1,), 2, 1)
            @test_throws ArgumentError run!(model, [100.0]; progress=false)                 # no coupling step
            @test_throws ArgumentError column(NACL_SYSTEM, (solution=1,), 2, 1; species=["Br-"])
            @test_throws ArgumentError column(NACL_SYSTEM, (solution=1,), 2, 1; species=[Species("Na+", 2, 1.0e-9)])
            @test_throws ArgumentError column(NACL_SYSTEM, Dict("Na+" => 1.0), 2, 1)        # numbers of definitions
            # Reservoirs given by a solution require geochemical reactions, and vice versa
            electrolyte = Dict("Na+" => 10.0, "Cl-" => 10.0)
            species = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9)]
            @test_throws ArgumentError Model1D(;
                length=0.01,
                grid,
                species,
                initial=electrolyte,
                left=Reservoir(; solution=2, potential=0.5),
                right=Reservoir(electrolyte; potential=0.0),
            )
            @test_throws ArgumentError Model1D(;
                length=0.01,
                grid,
                species,
                initial=electrolyte,
                left=Reservoir(electrolyte; potential=0.5),
                right=Reservoir(electrolyte; potential=0.0),
                database,
                input=NACL_SYSTEM,
            )
            transport = Model1D(;
                length=0.01,
                grid,
                species,
                initial=electrolyte,
                left=Reservoir(electrolyte; potential=0.5),
                right=Reservoir(electrolyte; potential=0.0),
            )
            @test_throws ArgumentError run!(transport, [100.0]; coupling_step=10.0, progress=false)
            @test_throws ArgumentError selected_output(transport, 0.0)
        end
    end
end
