# Aqueous species transported by the Nernst–Planck equations.

@doc raw"""
    Species(name, z, D)

Aqueous species with charge number `z` and diffusion coefficient `D` in free water
(m² s⁻¹). The effective diffusion coefficient in the porous medium is ``n \tau D`` (see
[`PorousMedium`](@ref)).

# Example
```julia
Species("Na+", 1, 1.33e-9)
```
"""
struct Species
    name::String
    z::Int
    D::Float64
    function Species(name::AbstractString, z::Integer, D::Real)
        isempty(name) && throw(ArgumentError("The species name must not be empty"))
        D >= 0 || throw(ArgumentError("The diffusion coefficient of $name must be non-negative, got $D"))
        return new(String(name), Int(z), Float64(D))
    end
end

Base.show(io::IO, s::Species) = print(io, "Species(\"", s.name, "\", z = ", s.z, ", D = ", s.D, " m²/s)")

"""
    concentration_vector(species, concentrations) -> Vector{Float64}

Concentrations (mol m⁻³) of `species`, in the same order, from a dictionary that maps species
names to concentrations. Species missing from the dictionary have zero concentration; names
that do not correspond to any species raise an error.
"""
function concentration_vector(species::AbstractVector{Species}, concentrations::AbstractDict)
    names = [s.name for s in species]
    for key in keys(concentrations)
        String(key) in names || throw(ArgumentError("Unknown species \"$key\"; the species are $(join(names, ", "))"))
    end
    c = [Float64(get(concentrations, s.name, 0.0)) for s in species]
    all(>=(0), c) || throw(ArgumentError("Concentrations must be non-negative"))
    return c
end

@doc raw"""
    charge_imbalance(species, c) -> Float64

Relative charge imbalance ``\left|\sum_i z_i c_i\right| / \sum_i |z_i| c_i`` of a solution
(zero for an electroneutral solution and for a solution without ions).
"""
function charge_imbalance(species::AbstractVector{Species}, c::AbstractVector)
    total = sum(abs(s.z) * ci for (s, ci) in zip(species, c))
    total == 0 && return 0.0
    return abs(sum(s.z * ci for (s, ci) in zip(species, c))) / total
end
