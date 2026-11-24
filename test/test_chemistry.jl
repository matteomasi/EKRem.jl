# Tests of the high-level chemistry layer (PhreeqcRM), on calcite in a closed system. The
# checks follow from mass and charge balances and from equilibrium, not from database values.

using Test
using EKRem

# Solution 2 introduces Na and Cl as components, so that HCl can be added to the cells
const CALCITE_SYSTEM = """
SOLUTION 1  Pure water
    pH 7 charge
SOLUTION 2  NaCl
    pH 7 charge
    Na 1
    Cl 1
EQUILIBRIUM_PHASES 1
    Calcite 0 10
KINETICS 1
    Calcite
        -m0 0.01
        -parms 50 0.6
        -tol 1e-12
SELECTED_OUTPUT 1
    -reset false
    -pH true
    -saturation_indices Calcite
USER_PUNCH 1
    -headings calcite_moles
    10 PUNCH KIN("Calcite")
END
"""

@testset "Chemistry" begin
    lib = phreeqcrm_library()
    if lib === nothing
        @info "No PhreeqcRM library set (see `set_phreeqcrm_library`): chemistry tests skipped"
        @test_skip false
    else
        database = get(ENV, "EKREM_PHREEQC_DATABASE", joinpath(dirname(dirname(lib)), "database", "phreeqc.dat"))
        closed_system(ncells, initial) = Chemistry(ncells; database, input=CALCITE_SYSTEM, initial)
        total(chemistry, component) =
            component_concentrations(chemistry)[:, findfirst(==(component), chemistry.components)]

        @testset "Equilibrium with calcite" begin
            chemistry = closed_system(2, (solution=1, equilibrium_phases=1))
            species = aqueous_species(chemistry)
            c = species_concentrations(chemistry)
            z = [s.z for s in species]
            @test [s.name for s in species] == chemistry.species
            @test size(c) == (2, length(species))
            @test all(c .>= 0)
            @test all(abs.(c * z) .< 1.0e-8 .* (c * abs.(z)))                 # electroneutrality
            @test all(abs.(selected_output(chemistry)["si_Calcite"]) .< 1.0e-8)
            @test total(chemistry, "Ca") ≈ total(chemistry, "C") rtol = 1.0e-8   # CaCO3 dissolved

            # Transferring the species concentrations back and equilibrating again leaves the
            # solution unchanged (species above 1 µmol/L; redox species are undetermined in
            # a solution without a redox couple)
            set_species_concentrations!(chemistry, c)
            react!(chemistry, 0.0)
            major = vec(any(c .> 1.0e-3; dims=1))
            @test species_concentrations(chemistry)[:, major] ≈ c[:, major] rtol = 1.0e-8
            close(chemistry)
            close(chemistry)   # closing twice is allowed
            @test chemistry.id == -1
        end

        @testset "Transfer of species concentrations" begin
            # HCl is added to the species concentrations, as transport would do
            amounts = [0.0, 1.0, 2.0, 5.0]   # mol/m³
            chemistry = closed_system(4, (solution=1, equilibrium_phases=1))
            ca0 = total(chemistry, "Ca")
            c = species_concentrations(chemistry)
            for name in ("H+", "Cl-")
                c[:, findfirst(==(name), chemistry.species)] .+= amounts
            end
            set_species_concentrations!(chemistry, c)
            react!(chemistry, 0.0)
            ca = total(chemistry, "Ca")
            @test total(chemistry, "Cl") ≈ amounts atol = 1.0e-12                # mass of Cl conserved
            @test ca ≈ total(chemistry, "C") rtol = 1.0e-8                         # CaCO3 dissolved
            @test all(abs.(selected_output(chemistry)["si_Calcite"]) .< 1.0e-8)
            # CaCO3 + 2 H+ → Ca+2 + CO2 + H2O and CaCO3 + H+ → Ca+2 + HCO3-: one mole of acid
            # dissolves between 1/2 and 1 mole of calcite
            dissolved = ca .- ca0
            @test all(amounts[2:end] / 2 .<= dissolved[2:end] .<= amounts[2:end])
            @test_throws ArgumentError set_species_concentrations!(chemistry, c[1:2, :])   # one row per cell
            close(chemistry)
        end

        @testset "Calcite kinetics" begin
            chemistry = closed_system(1, (solution=1, kinetics=1))
            ca = Float64[]
            for _ in 1:6
                react!(chemistry, 600.0)
                push!(ca, total(chemistry, "Ca")[1])
                # Dissolved calcite (mol per litre of water) appears as Ca in solution
                remaining = selected_output(chemistry)["calcite_moles"][1]
                @test 1000 * (0.01 - remaining) ≈ ca[end] rtol = 1.0e-8
            end
            @test chemistry.time == 3600.0
            @test issorted(ca) && ca[1] > 0
            @test selected_output(chemistry)["si_Calcite"][1] < 0             # still undersaturated
            # The result does not depend on the length of the reaction steps
            one_step = closed_system(1, (solution=1, kinetics=1))
            react!(one_step, 3600.0)
            @test total(one_step, "Ca")[1] ≈ ca[end] rtol = 1.0e-3
            @test_throws ArgumentError react!(one_step, -1.0)
        end

        @testset "Input checks" begin
            @test_throws ArgumentError closed_system(2, (solution=1, minerals=1))           # unknown entity
            @test_throws ArgumentError closed_system(2, (solution=[1, 1, 1],))              # one number per cell
            @test_throws ArgumentError Chemistry(1; database, input="missing.pqi", initial=(solution=1,))
            @test_throws ArgumentError Chemistry(1; database="missing.dat", input=CALCITE_SYSTEM, initial=(solution=1,))
            # Invalid PHREEQC input (PHREEQC also prints the error message on the console)
            invalid = "SOLUTION 1\n    pH seven\nEND\n"
            redirect_stdio(; stdout=devnull, stderr=devnull) do
                @test_throws PhreeqcRMError Chemistry(1; database, input=invalid, initial=(solution=1,))
            end
        end
    end
end
