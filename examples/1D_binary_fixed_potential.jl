#=
# Binary electrolyte under a fixed voltage

A 1 cm column of a water-saturated porous medium (porosity 0.4, tortuosity factor 0.5) lies
between two reservoirs with NaCl at 10 mol/m³ (anolyte) and 100 mol/m³ (catholyte). The
potential is fixed at 0.5 V at the anode and at 0 V at the cathode. The model is run until
the steady state is reached, and the result is compared with the analytical steady state

```math
c(x) = c_L + (c_R - c_L) \frac{x}{L}, \qquad
\phi(x) = \phi_L + (\phi_R - \phi_L) \frac{\ln(c(x) / c_L)}{\ln(c_R / c_L)}.
```

At steady state, the flux of each ion,
``J_i = -D_i \left( \mathrm{d}c_i/\mathrm{d}x + z_i c_i \, \mathrm{d}\psi/\mathrm{d}x \right)``
with ``\psi = F \phi / (R T)``, is constant along the column. Electroneutrality requires
``c_+ = c_- = c``, so that

```math
\frac{J_+}{D_+} = -\left( \frac{\mathrm{d}c}{\mathrm{d}x} + c \frac{\mathrm{d}\psi}{\mathrm{d}x} \right),
\qquad
\frac{J_-}{D_-} = -\left( \frac{\mathrm{d}c}{\mathrm{d}x} - c \frac{\mathrm{d}\psi}{\mathrm{d}x} \right).
```

The sum of the two equations shows that ``\mathrm{d}c/\mathrm{d}x`` is constant, so the
concentration is linear. Their difference shows that ``c \, \mathrm{d}\psi/\mathrm{d}x`` is
constant; with a linear ``c``, its integration gives the logarithmic potential, whose
constant is fixed by the boundary potentials. The diffusion coefficients enter only through
the constants ``J_i / D_i``, so the steady state depends neither on them nor on the porosity
and the tortuosity factor. The comparison is repeated on three grids to show the
second-order convergence of the discretisation.

The example script can be run from the repository root:

```console
julia --project=examples examples/1D_binary_fixed_potential.jl
```
=#

using EKRem
using Plots
using Printf
default(; guidefontsize=13, tickfontsize=11, legendfontsize=10)   # titles: default size

output = joinpath(@__DIR__, "output") #src
mkpath(output) #src

L = 0.01                   # column length (m)
cL, cR = 10.0, 100.0       # NaCl concentrations of the reservoirs (mol/m³)
φL, φR = 0.5, 0.0          # potentials of the reservoirs (V)
t_end = 1.0e6              # s; much longer than the diffusion time L² / (τ D) ≈ 1.5e5 s

electrolyte(c) = Dict("Na+" => c, "Cl-" => c)
species = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9)]

c_exact(x) = cL + (cR - cL) * x / L
φ_exact(x) = φL + (φR - φL) * log(c_exact(x) / cL) / log(cR / cL)
#md nothing #hide

# Simulations on three grids
grids = (25, 50, 100)                     # number of nodes
models = Model1D[]
println("Simulations")
for nodes in grids
    model = Model1D(;
        length=L,
        grid=collect(range(0, L; length=nodes)),
        species,
        initial=electrolyte(cL),
        left=Reservoir(electrolyte(cL); potential=φL),
        right=Reservoir(electrolyte(cR); potential=φR),
        porosity=0.4,
        tortuosity=0.5,
    )
    println(" grid with $nodes nodes:")
    run!(model, [t_end])
    push!(models, model)
end

# Errors with respect to the analytical solution (maximum over the nodes)
spacing = [L / (nodes - 1) for nodes in grids]
errors_c = [maximum(abs.(concentration(m, "Na+", t_end) .- c_exact.(coordinates(m)))) / cR for m in models]
errors_φ = [maximum(abs.(electric_potential(m, t_end) .- φ_exact.(coordinates(m)))) / (φL - φR) for m in models]

println("\nResults\nNodes   h (mm)   max error of c / c_R   max error of φ / Δφ")
for (nodes, h, err_c, err_φ) in zip(grids, spacing, errors_c, errors_φ)
    @printf("%5d   %6.3f   %20.2e   %20.2e\n", nodes, h * 1e3, err_c, err_φ)
end

# Order of convergence ``p`` with respect to the grid spacing, assuming error ``\propto h^p``,
# estimated from two successive grids (the spacing is not exactly halved, hence the general
# formula)
order(e1, e2, h1, h2) = log(e1 / e2) / log(h1 / h2)
println("\nOrder p of convergence with respect to the grid spacing h, assuming error ∝ hᵖ,")
println("from successive grids, p = ln(e₁/e₂) / ln(h₁/h₂):")
for k in 1:(length(spacing) - 1)
    @printf(
        "  h = %.3f → %.3f mm:  concentration p = %.2f, potential p = %.2f\n",
        spacing[k] * 1e3,
        spacing[k + 1] * 1e3,
        order(errors_c[k], errors_c[k + 1], spacing[k], spacing[k + 1]),
        order(errors_φ[k], errors_φ[k + 1], spacing[k], spacing[k + 1])
    )
end
println("(second-order scheme: p → 2 as h → 0)\n")

# Profiles on the finest grid
model = models[end]
x = coordinates(model)
xf = range(0, L; length=400)
marks = 1:5:length(x)
p1 = plot(; title="Concentration", xlabel="x (mm)", ylabel="c (mol/m³)", legend=:topleft)
scatter!(
    p1,
    x[marks] .* 1e3,
    concentration(model, "Na+", t_end)[marks];
    label="EKRem, Na⁺",
    marker=:square,
    markercolor=:white,
    markersize=6,
)
scatter!(
    p1,
    x[marks] .* 1e3,
    concentration(model, "Cl-", t_end)[marks];
    label="EKRem, Cl⁻",
    marker=:circle,
    markercolor=:black,
    markersize=6,
)
plot!(p1, xf .* 1e3, c_exact.(xf); color=:black, label="analytical")
p2 = plot(; title="Electric potential", xlabel="x (mm)", ylabel="φ (V)", legend=:topright)
scatter!(
    p2,
    x[marks] .* 1e3,
    electric_potential(model, t_end)[marks];
    label="EKRem",
    markercolor=:black,
    markerstrokewidth=0,
    markersize=6,
)
plot!(p2, xf .* 1e3, φ_exact.(xf); color=:black, label="analytical")
figure = plot(p1, p2; layout=(1, 2), size=(1000, 400), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_binary_fixed_potential.png")) #src
println("Figure saved in ", joinpath(output, "1D_binary_fixed_potential.png")) #src
#md figure
