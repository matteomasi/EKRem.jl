#=
EKRem example: benchmark with the Krabbenhøft analytical solution for the migration of
Na⁺, Cl⁻, K⁺ and OH⁻ under a fixed electric potential difference.

Setup

The domain is a 1D column of length L = 10 mm, filled with free solution (porosity and
tortuosity equal to 1). The electric potential is fixed at 3 V at the anode (x = 0) and at
0 V at the cathode (x = L), and is described by the Poisson equation. The concentrations at
the boundaries are (mol/m³):

    x = 0 (anode):     Na⁺ 25,   Cl⁻ 0,     K⁺ 83,   OH⁻ 108
    x = L (cathode):   Na⁺ 525,  Cl⁻ 500,   K⁺ 83,   OH⁻ 108

The column initially has the composition of the cathode boundary. Chloride migrates towards
the anode, and the other ions are displaced by the Coulombic interactions. The input data
are taken from Table 3 in Sprocati et al. (2019). The concentrations and the electric
potential are shown after 0.25 h and 1.5 h, and at steady state (24 h).

Analytical solution

At steady state, the electric potential is (Krabbenhøft and Krabbenhøft, 2008; Eq. 16 in
Sprocati et al., 2019):

    φ(x) = φ₀ + ln(1 + (Δc / c₀)(x / L)) / ln(c_L / c₀) · Δφ

and the concentration of each ion is (Eq. 17):

    cᵢ±(x) = (cᵢ,₀± c_L^(±α) − cᵢ,L± c₀^(±α)) / (c₀ c_L^(±α) − c_L c₀^(±α)) · ((x / L) Δc + c₀)
           − (cᵢ,₀± c_L − cᵢ,L± c₀) / (c₀ c_L^(±α) − c_L c₀^(±α)) · ((x / L) Δc + c₀)^(±α)

where:
- φ₀ and φ_L are the electric potentials at x = 0 and x = L, and Δφ = φ_L − φ₀;
- c₀ and c_L are the concentrations of all cations (or anions) at x = 0 and x = L,
  and Δc = c_L − c₀;
- cᵢ,₀± and cᵢ,L± are the concentrations of ion i at x = 0 and x = L;
- α = k Δφ / ln(c₀ / c_L), with k = F / (R T);
- the upper sign applies to cations and the lower sign to anions.


References

Krabbenhøft, K., Krabbenhøft, J., 2008. Application of the Poisson–Nernst–Planck equations
to the migration test. Cem. Concr. Res. 38, 77–88.
https://doi.org/10.1016/j.cemconres.2007.08.006

Sprocati, R., Masi, M., Muniruzzaman, M., Rolle, M., 2019. Modeling electrokinetic transport
and biogeochemical reactions in porous media: A multidimensional Nernst–Planck–Poisson
approach with PHREEQC coupling. Adv. Water Resour. 127, 134–147.
https://doi.org/10.1016/j.advwatres.2019.03.011

Run from the repository root:
    julia --project=examples examples/1D_benchmark_krabbenhoft.jl
=#

using EKRem
using Plots
using Printf

output = joinpath(@__DIR__, "output")
mkpath(output)

L = 0.01                   # length of the domain (m)
φ0, φL = 3.0, 0.0          # electric potential at x = 0 (anode) and at x = L (cathode) (V)
T = 298.15                 # temperature (K)

# Aqueous diffusion coefficients Dᵢᵃᑫ (m²/s) and boundary concentrations cᵢ,₀ and cᵢ,L (mol/m³)
species = [
    Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), Species("K+", 1, 1.96e-9), Species("OH-", -1, 5.27e-9)
]
ci0 = Dict("Na+" => 25.0, "Cl-" => 0.0, "K+" => 83.0, "OH-" => 108.0)
ciL = Dict("Na+" => 525.0, "Cl-" => 500.0, "K+" => 83.0, "OH-" => 108.0)

# Analytical steady state (Eqs. 16 and 17)
k = FARADAY / (GAS_CONSTANT * T)
c0 = sum(ci0[sp.name] for sp in species if sp.z > 0)   # concentration of all cations or anions at x = 0
cL = sum(ciL[sp.name] for sp in species if sp.z > 0)   # concentration of all cations or anions at x = L
Δc = cL - c0
Δφ = φL - φ0
α = k * Δφ / log(c0 / cL)

φ_exact(x) = φ0 + log(1 + Δc / c0 * x / L) / log(cL / c0) * Δφ

# Coefficients of Eq. (17), cᵢ = Aᵢ c + Bᵢ (c / m)^(±α), with the powers divided by m^(±α)
function coefficients(sp::Species)
    p = sp.z * α                  # ±α
    m = p > 0 ? cL : c0           # (c / m)^p ≤ 1 for c₀ ≤ c ≤ c_L
    a0, aL = ci0[sp.name], ciL[sp.name]
    den = c0 * (cL / m)^p - cL * (c0 / m)^p
    A = (a0 * (cL / m)^p - aL * (c0 / m)^p) / den
    B = -(a0 * cL - aL * c0) / den
    return A, B, p, m
end
function c_exact(sp::Species, x)
    A, B, p, m = coefficients(sp)
    c = x / L * Δc + c0
    return A * c + B * (c / m)^p
end
# Current density J = F Σ zᵢ Nᵢ, with the ionic flux Nᵢ = −Dᵢ (Δc / L) Aᵢ (1 − zᵢ α) at steady state
J_exact = FARADAY * sum(sp.z * (-sp.D * Δc / L * coefficients(sp)[1] * (1 - sp.z * α)) for sp in species)

model = Model1D(;
    length=L,
    grid=(h_min=1.25e-6, h_max=1.0e-5, growth=1.1),
    species,
    initial=ciL,
    left=Reservoir(ci0; potential=φ0),
    right=Reservoir(ciL; potential=φL),
    potential=Poisson(),
    temperature=T,
)
display(model)
println()

# Transient profiles after 0.25 h and 1.5 h; steady state at 24 h, longer than the diffusion
# time L² / Dᵢ ≈ 21 h of Na⁺
hour = 3600.0
transient = [0.25, 1.5] .* hour
t_end = 24 * hour
println("\nSimulation")
run!(model, [transient; t_end])

x = coordinates(model)
J_end = current_density(model, t_end)[1]   # between the first two nodes
errors_c = [maximum(abs.(concentration(model, sp.name, t_end) .- c_exact.(Ref(sp), x))) for sp in species]
error_φ = maximum(abs.(electric_potential(model, t_end) .- φ_exact.(x))) / abs(Δφ)

println("\nResults at t = $(t_end / hour) h (α = $(round(α; digits=2)))")
for (sp, err) in zip(species, errors_c)
    @printf("  max error of c(%s) / c_L: %.2e\n", sp.name, err / cL)
end
@printf("  max error of φ / |Δφ|: %.2e\n", error_φ)
@printf("  current density J: %.2f A/m² (analytical steady state %.2f A/m²)\n", J_end, J_exact)

# Concentration and electric potential at the transient times (dashed lines) and at steady
# state (markers, at the nodes closest to equally spaced points)
xf = sort(unique([range(0, L; length=2001); range(0, 2e-4; length=401); range(L - 1e-3, L; length=401)]))
marks = unique([argmin(abs.(x .- ξ)) for ξ in range(0, L; length=41)])
colors = palette(:viridis, length(transient) + 1)

function profile_panel(exact, numerical, title, ylabel; legend=false)
    p = plot(; title, xlabel="x (mm)", ylabel, legend)
    for (j, t) in enumerate(transient)
        plot!(p, x .* 1e3, numerical(t); color=colors[j], linestyle=:dash, label="EKRem, t = $(t / hour) h")
    end
    scatter!(
        p,
        x[marks] .* 1e3,
        numerical(t_end)[marks];
        label="EKRem, t = $(Int(t_end / hour)) h",
        markercolor=:black,
        markerstrokewidth=0,
        markersize=2.5,
    )
    plot!(p, xf .* 1e3, exact.(xf); color=:black, label="analytical (steady state)")
    return p
end

panels = [
    profile_panel(
        ξ -> c_exact(sp, ξ),
        t -> concentration(model, sp.name, t),
        sp.name,
        "C (mol/m³)";
        legend=sp == first(species) ? :topleft : false,
    ) for sp in species
]
push!(panels, profile_panel(φ_exact, t -> electric_potential(model, t), "Electric potential", "φ (V)"))

figure = plot(panels...; layout=@layout([a b c; d e _]), size=(1300, 750), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_benchmark_krabbenhoft.png"))
println("Figure saved in ", joinpath(output, "1D_benchmark_krabbenhoft.png"))
