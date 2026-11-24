# Tests for the low-level PhreeqcRM C interface.

using Test
using TOML
using EKRem
using EKRem.LibPhreeqcRM

const CALCITE_INPUT = """
SOLUTION 1  Pure water
    temp    25
    pH      7 charge
EQUILIBRIUM_PHASES 1
    Calcite 0 10
SELECTED_OUTPUT 1
    -reset              false
    -pH                 true
    -charge_balance     true
    -totals             Ca C
    -saturation_indices Calcite
END
"""

@testset "PhreeqcRM library preference" begin
    # The path is saved in LocalPreferences.toml next to the given Project.toml. A temporary
    # project is used, so that no real project is modified.
    mktempdir() do dir
        uuid = string(Base.PkgId(EKRem).uuid)
        project = joinpath(dir, "Project.toml")
        write(project, "[deps]\nEKRem = \"$uuid\"\n")
        library = joinpath(dir, "lib", "PhreeqcRM.dll")
        LibPhreeqcRM._save_library(library; project_toml=project)
        prefs = TOML.parsefile(joinpath(dir, "LocalPreferences.toml"))
        @test prefs["EKRem"]["phreeqcrm_library"] == library
        # Saving again replaces the previous value
        LibPhreeqcRM._save_library("other.dll"; project_toml=project)
        @test TOML.parsefile(joinpath(dir, "LocalPreferences.toml"))["EKRem"]["phreeqcrm_library"] == "other.dll"
        # The project file itself is not modified
        @test read(project, String) == "[deps]\nEKRem = \"$uuid\"\n"
    end
end

@testset "PhreeqcRM simulation: calcite in pure water" begin
    lib = phreeqcrm_library()
    if lib === nothing
        @info "No PhreeqcRM library set (see `set_phreeqcrm_library`): PhreeqcRM tests skipped"
        @test_skip false
    else
        database = get(ENV, "EKREM_PHREEQC_DATABASE", joinpath(dirname(dirname(lib)), "database", "phreeqc.dat"))
        @test isfile(database)

        mktempdir() do dir
            nxyz = 5
            id = LibPhreeqcRM.RM_Create(nxyz, 1)
            @test id >= 0
            try
                LibPhreeqcRM.RM_SetFilePrefix(id, joinpath(dir, "calcite"))
                LibPhreeqcRM.RM_OpenFiles(id)
                @test LibPhreeqcRM.RM_GetGridCellCount(id) == nxyz
                LibPhreeqcRM.RM_LoadDatabase(id, database)
                LibPhreeqcRM.RM_SetSelectedOutputOn(id, true)
                LibPhreeqcRM.RM_SetSpeciesSaveOn(id, true)   # before RM_FindComponents

                # Definitions in the worker instances (SELECTED_OUTPUT) and in InitialPhreeqc
                input = joinpath(dir, "calcite.pqi")
                write(input, CALCITE_INPUT)
                LibPhreeqcRM.RM_RunFile(id, true, true, false, input)
                ncomps = LibPhreeqcRM.RM_FindComponents(id)
                comps = [LibPhreeqcRM.RM_GetComponent(id, i) for i in 0:(ncomps - 1)]
                @test ["H", "O", "Charge", "C", "Ca"] ⊆ comps

                # Solution 1 and equilibrium phases 1 in every cell; reaction step
                ic1 = fill(-1, nxyz, 7)
                ic1[:, 1] .= 1
                ic1[:, 2] .= 1
                LibPhreeqcRM.RM_InitialPhreeqc2Module(id, ic1)
                LibPhreeqcRM.RM_RunCells(id)

                LibPhreeqcRM.RM_SetCurrentSelectedOutputUserNumber(id, 1)
                so = LibPhreeqcRM.RM_GetSelectedOutput(id)
                headings = [LibPhreeqcRM.RM_GetSelectedOutputHeading(id, j) for j in 0:(size(so, 2) - 1)]
                @test size(so, 1) == nxyz
                column(name) = so[:, findfirst(==(name), headings)]
                pH = column("pH")
                charge = column("charge(eq)")
                ca = column("Ca(mol/kgw)")
                c = column("C(mol/kgw)")
                si = column("si_Calcite")
                @info "Calcite in pure water" pH = pH[1] Ca = ca[1] C = c[1] si_Calcite = si[1]

                # Identical cells give identical results
                @test all(so .== so[1:1, :])
                # Equilibrium with calcite
                @test all(abs.(si) .< 1.0e-8)
                # Dissolution of CaCO3 in a closed system: equal amounts of Ca and C
                @test all(ca .> 0)
                @test ca ≈ c rtol = 1.0e-10
                # Electroneutrality of the equilibrated solution
                @test all(abs.(charge) .< 1.0e-8 .* 2ca)
                # Dissolution of a carbonate mineral raises the pH above neutrality
                @test all(pH .> 7)

                @testset "Aqueous species" begin
                    nspecies = LibPhreeqcRM.RM_GetSpeciesCount(id)
                    species = [LibPhreeqcRM.RM_GetSpeciesName(id, i) for i in 0:(nspecies - 1)]
                    z = LibPhreeqcRM.RM_GetSpeciesZ(id)
                    d25 = LibPhreeqcRM.RM_GetSpeciesD25(id)
                    conc = LibPhreeqcRM.RM_GetSpeciesConcentrations(id)   # mol/L
                    @test length(z) == length(d25) == nspecies
                    @test size(conc) == (nxyz, nspecies)
                    @test allunique(species)

                    # The main species of the Ca–C–H–O system are present
                    @test ["H+", "OH-", "Ca+2", "HCO3-", "CO3-2"] ⊆ species

                    # The charge numbers agree with the charge written in the species names
                    # (e.g. "CO3-2" → −2, "CaHCO3+" → +1, "CaCO3" → 0)
                    function charge_from_name(name)
                        m = match(r"([+-])(\d*)$", name)
                        m === nothing && return 0
                        return (m[1] == "+" ? 1 : -1) * (isempty(m[2]) ? 1 : parse(Int, m[2]))
                    end
                    @test z == charge_from_name.(species)

                    # Diffusion coefficients are non-negative; H+ has the largest one
                    @test all(d25 .>= 0)
                    @test species[argmax(d25)] == "H+"

                    # Concentrations are non-negative, and the solution is electroneutral:
                    # Σ zᵢ cᵢ = 0 relative to Σ |zᵢ| cᵢ
                    @test all(conc .>= 0)
                    @test all(abs.(conc * z) .< 1.0e-10 .* (conc * abs.(z)))

                    # The H+ concentration is consistent with the pH of the selected output
                    # (they differ only by the activity coefficient, close to 1 in this dilute
                    # solution)
                    h = conc[:, findfirst(==("H+"), species)]
                    @test all(abs.(-log10.(h) .- pH) .< 0.05)
                end

                LibPhreeqcRM.RM_CloseFiles(id)
            finally
                LibPhreeqcRM.RM_Destroy(id)
            end
        end
    end
end
