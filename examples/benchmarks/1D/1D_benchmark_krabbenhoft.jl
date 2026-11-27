#=
# Krabbenhøft and Krabbenhøft benchmark

Electromigration of Na⁺, Cl⁻, K⁺ and OH⁻ under a fixed electric potential difference,
compared with the analytical steady state of Krabbenhøft and Krabbenhøft (2008).

## Setup

The domain is a column of length ``L`` = 10 mm, filled with free solution (porosity and
tortuosity factor equal to 1). The electric potential is fixed at 3 V at the anode
(``x = 0``) and at 0 V at the cathode (``x = L``), and it is described by the Poisson
equation. The concentrations at the boundaries are (mol/m³):

| Boundary            | Na⁺ | Cl⁻ | K⁺ | OH⁻ |
|:--------------------|----:|----:|---:|----:|
| ``x = 0`` (anode)   |  25 |   0 | 83 | 108 |
| ``x = L`` (cathode) | 525 | 500 | 83 | 108 |

The column initially has the composition of the cathode boundary. Chloride migrates towards
the anode, and the other ions are displaced by the Coulombic interactions. The input data
are taken from Table 3 in Sprocati et al. (2019). The concentrations and the electric
potential are shown after 0.25 h and 1.5 h, and at steady state (24 h).

## Analytical solution

At steady state, the electric potential is (Krabbenhøft and Krabbenhøft, 2008; Eq. 16 in
Sprocati et al., 2019)

```math
\phi(x) = \phi_0 + \frac{\ln\left(1 + \dfrac{\Delta c}{c_0} \dfrac{x}{L}\right)}{\ln(c_L / c_0)} \, \Delta\phi,
```

and the concentration of each ion is (Eq. 17)

```math
c_i^\pm(x) = \frac{c_{i,0}^\pm c_L^{\pm\alpha} - c_{i,L}^\pm c_0^{\pm\alpha}}{c_0 c_L^{\pm\alpha} - c_L c_0^{\pm\alpha}}
\left( \frac{x}{L} \Delta c + c_0 \right)
- \frac{c_{i,0}^\pm c_L - c_{i,L}^\pm c_0}{c_0 c_L^{\pm\alpha} - c_L c_0^{\pm\alpha}}
\left( \frac{x}{L} \Delta c + c_0 \right)^{\pm\alpha},
```

where:
- ``\phi_0`` and ``\phi_L`` are the electric potentials at ``x = 0`` and ``x = L``, and
  ``\Delta\phi = \phi_L - \phi_0``;
- ``c_0`` and ``c_L`` are the concentrations of all cations (or anions) at ``x = 0`` and
  ``x = L``, and ``\Delta c = c_L - c_0``;
- ``c_{i,0}^\pm`` and ``c_{i,L}^\pm`` are the concentrations of ion ``i`` at ``x = 0`` and
  ``x = L``;
- ``\alpha = k \Delta\phi / \ln(c_0 / c_L)``, with ``k = F / (R T)``;
- the upper sign applies to cations and the lower sign to anions.


The example script can be run from the repository root:

```console
julia --project=examples examples/benchmarks/1D/1D_benchmark_krabbenhoft.jl
```
=#

using SimEK
using Plots
using Printf
default(; guidefontsize=13, tickfontsize=11, legendfontsize=10)   # titles: default size

output = joinpath(@__DIR__, "output") #src
mkpath(output) #src

L = 0.01                   # length of the domain (m)
φ0, φL = 3.0, 0.0          # electric potential at x = 0 (anode) and at x = L (cathode) (V)
T = 298.15                 # temperature (K)
#md nothing #hide

# Aqueous diffusion coefficients ``D_i`` (m²/s) and boundary concentrations ``c_{i,0}`` and
# ``c_{i,L}`` (mol/m³)
species = [
    Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), Species("K+", 1, 1.96e-9), Species("OH-", -1, 5.27e-9)
]
ci0 = Dict("Na+" => 25.0, "Cl-" => 0.0, "K+" => 83.0, "OH-" => 108.0)
ciL = Dict("Na+" => 525.0, "Cl-" => 500.0, "K+" => 83.0, "OH-" => 108.0)
#md nothing #hide

# Analytical steady state (Eqs. 16 and 17)
k = FARADAY / (GAS_CONSTANT * T)
c0 = sum(ci0[sp.name] for sp in species if sp.z > 0)   # concentration of all cations or anions at x = 0
cL = sum(ciL[sp.name] for sp in species if sp.z > 0)   # concentration of all cations or anions at x = L
Δc = cL - c0
Δφ = φL - φ0
α = k * Δφ / log(c0 / cL)

φ_exact(x) = φ0 + log(1 + Δc / c0 * x / L) / log(cL / c0) * Δφ
#md nothing #hide

# Coefficients of Eq. (17), ``c_i = A_i c + B_i (c / m)^{\pm\alpha}``, with the powers divided
# by ``m^{\pm\alpha}``
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
#md nothing #hide
# Current density ``J = F \sum_i z_i N_i``, with the ionic flux
# ``N_i = -D_i (\Delta c / L) A_i (1 - z_i \alpha)`` at steady state
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
display(model) #src
println() #src
#md model

# Transient profiles after 0.25 h and 1.5 h; steady state at 24 h, longer than the diffusion
# time ``L^2 / D_i \approx`` 21 h of Na⁺
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

# Concentration and electric potential at the transient times and at steady state (markers,
# at the nodes closest to equally spaced points)
xf = sort(unique([range(0, L; length=2001); range(0, 2e-4; length=401); range(L - 1e-3, L; length=401)]))
marks = unique([argmin(abs.(x .- ξ)) for ξ in range(0, L; length=41)])
colors = ["#0072B2", "#D55E00", "#009E73"]   # Okabe–Ito palette (colour-blind safe)

function profile_panel(exact, numerical, title, ylabel; legend=false)
    p = plot(; title, xlabel="x (mm)", ylabel, legend)
    for (j, t) in enumerate([transient; t_end])
        scatter!(
            p,
            x[marks] .* 1e3,
            numerical(t)[marks];
            label="SimEK, t = $(t / hour) h",
            markercolor=t == t_end ? :black : colors[j],
            markerstrokewidth=0,
            markersize=3.5,
        )
    end
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

figure = plot(panels...; layout=@layout([a b; c d; e _]), size=(1000, 1200), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_benchmark_krabbenhoft.png")) #src
println("Figure saved in ", joinpath(output, "1D_benchmark_krabbenhoft.png")) #src
#md figure

#=
## References

Krabbenhøft, K., Krabbenhøft, J., 2008. Application of the Poisson–Nernst–Planck equations
to the migration test. Cem. Concr. Res. 38, 77–88.
[https://doi.org/10.1016/j.cemconres.2007.08.006](https://doi.org/10.1016/j.cemconres.2007.08.006)

Sprocati, R., Masi, M., Muniruzzaman, M., Rolle, M., 2019. Modeling electrokinetic transport
and biogeochemical reactions in porous media: A multidimensional Nernst–Planck–Poisson
approach with PHREEQC coupling. Adv. Water Resour. 127, 134–147.
[https://doi.org/10.1016/j.advwatres.2019.03.011](https://doi.org/10.1016/j.advwatres.2019.03.011)
=#
