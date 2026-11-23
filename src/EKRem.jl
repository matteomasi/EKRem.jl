"""
    EKRem

Simulation of electrokinetic remediation. Multicomponent transport is described by the
Nernst–Planck equations, the electric potential by the Poisson or the charge-conservation
equation, and geochemical reactions are computed by PhreeqcRM.

Models are organised by spatial dimension (`Model1D`, `Model2D`, `Model3D`).
"""
module EKRem

using VoronoiFVM: fbernoulli

include("constants.jl")
include("physics/fluxes.jl")
include("phreeqcrm/capi.jl")

using .LibPhreeqcRM: set_phreeqcrm_library, phreeqcrm_library, PhreeqcRMError

export FARADAY, GAS_CONSTANT, VACUUM_PERMITTIVITY, WATER_RELATIVE_PERMITTIVITY
export debye_length, charge_relaxation_time
export sg_flux
export set_phreeqcrm_library, phreeqcrm_library, PhreeqcRMError

end # module EKRem
