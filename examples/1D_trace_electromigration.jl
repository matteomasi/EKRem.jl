#=
# Electromigration and electroosmosis of a trace ion

A 30 cm column contains a water-saturated porous medium (porosity ``n`` = 0.5, tortuosity
factor ``\tau`` = 0.7, electroosmotic permeability ``k_{eo}`` = 5 × 10⁻⁹ m² V⁻¹ s⁻¹). The
initial pore water contains NaCl at 100 mol/m³. The voltage applied at the anode is 6 V,
which corresponds to the gradient ``-\partial\phi/\partial x`` = 20 V/m. A trace cation T⁺
(0.01 mol/m³, ``D`` = 1.96 × 10⁻⁹ m²/s) enters from the left reservoir.

The electric field can be assumed to remain spatially uniform. Under this assumption, the
transport of the trace ion by diffusion, electromigration and electroosmosis is described
by the one-dimensional advection–diffusion equation

```math
n \frac{\partial c}{\partial t} = D^* \frac{\partial^2 c}{\partial x^2} - k_v \frac{\partial c}{\partial x},
```

where ``D^* = n \tau D`` is the effective diffusion coefficient and ``k_v`` accounts for
both the electromigration and the electroosmotic fluxes:

```math
k_v = -\frac{z F}{R T} D^* \left( \frac{\partial \phi}{\partial x} \right)
- k_{eo} \left( \frac{\partial \phi}{\partial x} \right).
```

With ``c(0, t) = c_{in}`` and ``c(x, 0) = 0``, the analytical solution on a semi-infinite
domain is (Ogata and Banks, 1961; Acar et al., 1990)

```math
c(x, t) = \frac{c_{in}}{2} \left[ \operatorname{erfc}\left( \frac{x - k_v t}{2 \sqrt{D^* t}} \right)
+ \exp\left( \frac{k_v x}{D^*} \right) \operatorname{erfc}\left( \frac{x + k_v t}{2 \sqrt{D^* t}} \right) \right],
```

evaluated at the time ``t/n``, because ``D^*`` and ``k_v`` refer to the bulk medium, whereas
the solute is stored in the pore water.

The column is long enough for the right boundary to have no influence at the output times.

## References

Acar, Y.B., Gale, R.J., Putnam, G.A., Hamed, J., Wong, R.L., 1990. Electrochemical processing
of soils: theory of pH gradient development by diffusion, migration, and linear convection.
J. Environ. Sci. Health A 25(6), 687–714. https://doi.org/10.1080/10934529009375590

Ogata, A., Banks, R.B., 1961. A solution of the differential equation of longitudinal
dispersion in porous media. U.S. Geological Survey Professional Paper 411-A, A1–A7.
https://doi.org/10.3133/pp411A

Run from the repository root:

```console
julia --project=examples examples/1D_trace_electromigration.jl
```
=#

using EKRem
using Plots
using Printf
using SpecialFunctions: erfc
default(; guidefontsize=13, tickfontsize=11, legendfontsize=10)   # titles: default size

output = joinpath(@__DIR__, "output") #src
mkpath(output) #src

L = 0.3                    # column length (m)
ΔV = 6.0                   # applied voltage (V)
cb, c_in = 100.0, 0.01     # background NaCl and inflow trace concentration (mol/m³)
porosity, tortuosity = 0.5, 0.7
k_eo = 5.0e-9              # electroosmotic permeability (m² V⁻¹ s⁻¹)
times = [10.0, 20.0, 30.0] .* 3600   # s

trace = Species("T+", 1, 1.96e-9)
model = Model1D(;
    length=L,
    grid=collect(range(0, L; length=601)),       # uniform, 0.5 mm
    species=[Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), trace],
    initial=Dict("Na+" => cb, "Cl-" => cb),
    left=Reservoir(Dict("Na+" => cb, "Cl-" => cb + c_in, "T+" => c_in); potential=ΔV),
    right=Reservoir(Dict("Na+" => cb, "Cl-" => cb); potential=0.0),
    porosity,
    tortuosity,
    electroosmotic_permeability=k_eo,
)
display(model) #src
println() #src
#md model

# The model is run up to the last output time.
elapsed = @elapsed run!(model, times)
@printf("Run time: %.1f s\n\n", elapsed)

Dstar = porosity * tortuosity * trace.D                     # effective diffusion coefficient D* (m²/s)
dφdx = -ΔV / L                                              # potential gradient (V/m)
kv = -trace.z * FARADAY / (GAS_CONSTANT * 298.15) * Dstar * dφdx - k_eo * dφdx   # k_v (m/s)
q = electroosmotic_flow(model, times[end])
@printf("Electroosmotic water flux q = %.3g m/s (pore-water velocity q/n = %.3g m/s)\n\n", q, q / porosity)
function exact(x, t)
    return c_in / 2 *
           (erfc((x - kv * t) / (2sqrt(Dstar * t))) + exp(kv * x / Dstar) * erfc((x + kv * t) / (2sqrt(Dstar * t))))
end
#md nothing #hide

# Position where the concentration falls to cᵢₙ/2 (linear interpolation between nodes)
function half_point(x, c)
    i = findlast(>=(c_in / 2), c)
    return x[i] + (c_in / 2 - c[i]) * (x[i + 1] - x[i]) / (c[i + 1] - c[i])
end

x = coordinates(model)
colors = ["#0072B2", "#D55E00", "#000000"]   # Okabe–Ito palette (colour-blind safe)
p1 = plot(; title="Trace ion T⁺", xlabel="x (cm)", ylabel="c/cᵢₙ", xlims=(-1, 31), legend=:topright)
xf = range(0, L; length=400)   # analytical solutions, drawn behind the points
for (k, t) in enumerate(times)
    plot!(p1, xf .* 100, exact.(xf, t / porosity) ./ c_in; color=:black, label=k == 1 ? "analytical" : "")
end
println("    t (s)   x(c = cᵢₙ/2), EKRem   x(c = cᵢₙ/2), analytical   max |c − c_exact| / cᵢₙ")
for (k, t) in enumerate(times)
    c = concentration(model, "T+", t)
    @printf(
        "%9.0f   %16.5f m   %21.5f m   %22.2e\n",
        t,
        half_point(x, c),
        half_point(x, exact.(x, t / porosity)),
        maximum(abs.(c .- exact.(x, t / porosity))) / c_in
    )
    scatter!(
        p1,
        x[1:10:end] .* 100,
        c[1:10:end] ./ c_in;
        label=@sprintf("EKRem, t = %.0f h", t / 3600),
        markercolor=colors[k],
        markerstrokewidth=0,
        markersize=6,
    )
end
# Electric potential at the last output time, compared with the uniform field assumed by the
# analytical solution
p2 = plot(; title="Electric potential", xlabel="x (cm)", ylabel="φ (V)", xlims=(-1, 31), legend=:topright)
plot!(p2, [0, L * 100], [ΔV, 0]; color=:black, label="uniform field")
scatter!(
    p2,
    x[1:30:end] .* 100,
    electric_potential(model, times[end])[1:30:end];
    label=@sprintf("EKRem, t = %.0f h", times[end] / 3600),
    markercolor=colors[end],
    markerstrokewidth=0,
    markersize=6,
)
figure = plot(p1, p2; layout=(1, 2), size=(1000, 400), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_trace_electromigration.png")) #src
println("\nFigure saved in ", joinpath(output, "1D_trace_electromigration.png")) #src
#md figure
