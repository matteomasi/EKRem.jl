"""
    EKRem

Simulation of electrokinetic remediation. Multicomponent transport is described by the
Nernst–Planck equations, the electric potential by the Poisson or the charge-conservation
equation, and geochemical reactions are computed by PhreeqcRM ([`Chemistry`](@ref)).

Models are organised by spatial dimension (`Model1D`, `Model2D`, `Model3D`).
"""
module EKRem

using ExtendableGrids: ExtendableGrids, ExtendableGrid, simplexgrid
using LinearSolve: UMFPACKFactorization
using OrdinaryDiffEqBDF: FBDF
using OrdinaryDiffEqRosenbrock: Rodas5P
using PrecompileTools: @compile_workload, @setup_workload
using Printf: @sprintf
using SciMLBase: SciMLBase
using SparseArrays: sparse
using VoronoiFVM: VoronoiFVM, fbernoulli

include("constants.jl")
include("electrolyte.jl")
include("species.jl")
include("porousmedium.jl")
include("grids.jl")
include("physics/fluxes.jl")
include("physics/potential.jl")
include("physics/boundaries.jl")
include("transport.jl")
include("phreeqcrm/capi.jl")
include("chemistry.jl")
include("coupling.jl")
include("models/model1d.jl")

using .LibPhreeqcRM: set_phreeqcrm_library, phreeqcrm_library, PhreeqcRMError

export FARADAY, GAS_CONSTANT, VACUUM_PERMITTIVITY, WATER_RELATIVE_PERMITTIVITY
export debye_length, conductivity, charge_relaxation_time
export Species, PorousMedium, grid_refined
export sg_flux
export PotentialModel, Poisson, ChargeConservation, FixedField
export BoundaryCondition, Reservoir, Compartment
export Model1D, run!, coordinates, output_times, concentration, electric_potential, voltage, current_density, amount
export electroosmotic_flow
export debye_ratio
export set_phreeqcrm_library, phreeqcrm_library, PhreeqcRMError
export Chemistry, aqueous_species, species_concentrations, set_species_concentrations!, component_concentrations
export react!, selected_output

# Precompilation workload: a quick simulation is run when the package is precompiled,
# to speed up the first simulation of a session.
@setup_workload begin
    function workload_model(n; potential=Poisson(), current=false)
        species = [Species("S$i", isodd(i) ? 1 : -1, 1.0e-9) for i in 1:n]
        # Equal concentrations of cations and anions (the last species is absent for odd n)
        c(value) = Dict(s.name => (iseven(n) || i < n ? value : 0.0) for (i, s) in enumerate(species))
        return Model1D(;
            length=1.0e-3,
            grid=(h_min=1.0e-4, h_max=2.0e-4),
            species=species,
            initial=c(1.0),
            left=current ? Reservoir(c(1.0); current_density=1.0) : Reservoir(c(1.0); potential=0.1),
            right=Reservoir(c(2.0); potential=0.0),
            potential,
        )
    end
    @compile_workload begin
        for n in 2:5
            run!(workload_model(n), [1.0]; progress=false)
        end
        model = workload_model(2)
        redirect_stdout(devnull) do
            run!(model, [1.0]; progress=true)   # the progress callback changes the compiled code
        end
        run!(model, [2.0]; solver=:FBDF, progress=false)
        concentration(model, "S1", 2.0)
        electric_potential(model, 2.0)
        voltage(model, 2.0)
        current_density(model, 2.0)
        amount(model, "S1", 2.0)
        sprint(show, MIME("text/plain"), model)
        run!(workload_model(2; potential=ChargeConservation(), current=true), [1.0]; progress=false)
        run!(workload_model(2; potential=FixedField()), [1.0]; progress=false)
        # Electrode compartment (water electrolysis)
        water = [Species("H+", 1, 1.0e-9), Species("OH-", -1, 1.0e-9)]
        c = Dict("H+" => 1.0, "OH-" => 1.0)
        compartment = Model1D(;
            length=1.0e-3,
            grid=(h_min=1.0e-4, h_max=2.0e-4),
            species=water,
            initial=c,
            left=Compartment(c; volume=1.0e-6, current_density=1.0),
            right=Reservoir(c; potential=0.0),
            potential=ChargeConservation(),
            area=1.0e-4,
        )
        run!(compartment, [1.0]; progress=false)
        # Compartments under voltage control
        compartments = Model1D(;
            length=1.0e-3,
            grid=(h_min=1.0e-4, h_max=2.0e-4),
            species=water,
            initial=c,
            left=Compartment(c; volume=1.0e-6, potential=0.1),
            right=Compartment(c; volume=1.0e-6, potential=0.0),
            potential=ChargeConservation(),
            area=1.0e-4,
        )
        run!(compartments, [1.0]; progress=false)
    end
end

end # module EKRem
