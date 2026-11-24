using Test

@testset "EKRem" begin
    include("test_constants.jl")
    include("test_electrolyte.jl")
    include("test_fluxes.jl")
    include("test_grids.jl")
    include("test_model1d.jl")
    include("test_phreeqcrm_capi.jl")
    include("test_chemistry.jl")
    include("test_coupling.jl")
    include("test_compartments.jl")
end
