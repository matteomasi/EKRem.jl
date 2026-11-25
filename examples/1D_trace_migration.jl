#=
# Migration of a trace ion

A 10 cm column contains a water-saturated porous medium (porosity 0.4, tortuosity factor
0.5). The initial pore water contains NaCl at 100 mol/m³. The voltage applied at the anode
is 1 V. A trace cation T⁺ (0.01 mol/m³) enters from the left reservoir.

The electric field can be assumed to remain spatially uniform. Under this assumption, the
transport of the trace ion is governed by the advection–dispersion equation

```math
\frac{\partial c}{\partial t} = \tau D \frac{\partial^2 c}{\partial x^2} - v \frac{\partial c}{\partial x},
\qquad v = \frac{\tau D z F E}{R T},
```

with ``c(0, t) = c_0`` and ``c(x, 0) = 0``. On a semi-infinite domain, its analytical
solution is (Ogata and Banks, 1961)

```math
\frac{c}{c_0} = \frac{1}{2} \left[ \operatorname{erfc}\left( \frac{x - v t}{2 \sqrt{\tau D t}} \right)
+ \exp\left( \frac{v x}{\tau D} \right) \operatorname{erfc}\left( \frac{x + v t}{2 \sqrt{\tau D t}} \right) \right].
```

The column is long enough for the right boundary to have no influence at the output times.

Run from the repository root:

```console
julia --project=examples examples/1D_trace_migration.jl
```
=#

using EKRem
using Plots
using Printf
using SpecialFunctions: erfc
default(; guidefontsize=13, tickfontsize=11, legendfontsize=10)   # titles: default size

output = joinpath(@__DIR__, "output") #src
mkpath(output) #src

L = 0.1                    # column length (m)
ΔV = 1.0                   # applied voltage (V)
cb, c0 = 100.0, 0.01       # background NaCl and trace concentration (mol/m³)
porosity, tortuosity = 0.4, 0.5
times = [2.5e4, 5.0e4, 1.0e5]   # s

trace = Species("T+", 1, 1.0e-9)
model = Model1D(;
    length=L,
    grid=collect(range(0, L; length=201)),       # uniform, 0.5 mm
    species=[Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), trace],
    initial=Dict("Na+" => cb, "Cl-" => cb),
    left=Reservoir(Dict("Na+" => cb, "Cl-" => cb + c0, "T+" => c0); potential=ΔV),
    right=Reservoir(Dict("Na+" => cb, "Cl-" => cb); potential=0.0),
    porosity,
    tortuosity,
)
display(model) #src
println() #src
#md model

# The model is run up to the last output time.
elapsed = @elapsed run!(model, times)
@printf("Run time: %.1f s\n\n", elapsed)

Dp = tortuosity * trace.D                                   # pore-water diffusion coefficient
v = Dp * FARADAY / (GAS_CONSTANT * 298.15) * trace.z * ΔV / L   # migration velocity (m/s)
exact(x, t) = c0 / 2 * (erfc((x - v * t) / (2sqrt(Dp * t))) + exp(v * x / Dp) * erfc((x + v * t) / (2sqrt(Dp * t))))
#md nothing #hide

# Position where the concentration falls to c₀/2 (linear interpolation between nodes)
function half_point(x, c)
    i = findlast(>=(c0 / 2), c)
    return x[i] + (c0 / 2 - c[i]) * (x[i + 1] - x[i]) / (c[i + 1] - c[i])
end

x = coordinates(model)
colors = ["#0072B2", "#D55E00", "#000000"]   # Okabe–Ito palette (colour-blind safe)
p = plot(;
    title="Trace ion T⁺",
    xlabel="x (cm)",
    ylabel="c/c₀",
    xlims=(0, 4),
    legend=:topright,
    size=(1000, 450),
    margin=5Plots.mm,
)
println("    t (s)   x(c = c₀/2), EKRem   x(c = c₀/2), analytical   max |c − c_exact| / c₀")
for (k, t) in enumerate(times)
    c = concentration(model, "T+", t)
    @printf(
        "%9.0f   %16.5f m   %21.5f m   %22.2e\n",
        t,
        half_point(x, c),
        half_point(x, exact.(x, t)),
        maximum(abs.(c .- exact.(x, t))) / c0
    )
    scatter!(
        p,
        x[1:2:end] .* 100,
        c[1:2:end] ./ c0;
        label=@sprintf("EKRem, t = %.0f h", t / 3600),
        markercolor=colors[k],
        markerstrokewidth=0,
        markersize=6,
    )
end
# Analytical solutions in front of the points
xf = range(0, 0.04; length=400)
for (k, t) in enumerate(times)
    plot!(p, xf .* 100, exact.(xf, t) ./ c0; color=:black, label=k == 1 ? "analytical" : "")
end
savefig(p, joinpath(output, "1D_trace_migration.png")) #src
println("\nFigure saved in ", joinpath(output, "1D_trace_migration.png")) #src
#md p
