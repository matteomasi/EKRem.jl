#=
EKRem example: migration of a trace ion in a uniform electric field

A 10 cm column contains water-saturated porous medium (porosity 0.4, tortuosity factor 0.5).
The initial pore water contains NaCl at 100 mol/m³. The voltage applied at the anode is 1 V. 
A trace cation T⁺ (0.01 mol/m³) enters from the left reservoir.
The electric field can be assumed to remain spatially uniform. Under this assumption, the transport 
of the trace ion is governed by the advection–dispersion equation:

    ∂c/∂t = τ D ∂²c/∂x² − v ∂c/∂x,    v = τ D z F E / (R T),

with c(0, t) = c₀ and c(x, 0) = 0. 
On a semi-infinite domain its analytical solution is the following (Ogata and Banks, 1961):

    c / c₀ = 0.5 * [erfc((x − v t) / (2 √(τ D t))) + exp(v x / (τ D)) erfc((x + v t) / (2 √(τ D t)))]

The column is long enough for the right boundary to have no influence at the considered output times.

Run from the repository root:
    julia --project=examples examples/1D_trace_migration.jl
=#

using EKRem
using Plots
using Printf
using SpecialFunctions: erfc

output = joinpath(@__DIR__, "output")
mkpath(output)

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
display(model)
println()
elapsed = @elapsed run!(model, times)
@printf("Run time: %.1f s\n\n", elapsed)

Dp = tortuosity * trace.D                                   # pore-water diffusion coefficient
v = Dp * FARADAY / (GAS_CONSTANT * 298.15) * trace.z * ΔV / L   # migration velocity (m/s)
exact(x, t) = c0 / 2 * (erfc((x - v * t) / (2sqrt(Dp * t))) + exp(v * x / Dp) * erfc((x + v * t) / (2sqrt(Dp * t))))

# Position where the concentration falls to c₀/2 (linear interpolation between nodes)
function half_point(x, c)
    i = findlast(>=(c0 / 2), c)
    return x[i] + (c0 / 2 - c[i]) * (x[i + 1] - x[i]) / (c[i + 1] - c[i])
end

x = coordinates(model)
colours = ["#0072B2", "#D55E00", "#009E73"]
markers = [:circle, :square, :diamond]
sizes = [5, 4.3, 6.0]
p = plot(;
    xlabel="x (cm)", ylabel="T⁺ concentration c/c₀", xlims=(0, 4), legend=:topright, size=(600, 400), margin=5Plots.mm
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
    plot!(p, x .* 100, exact.(x, t) ./ c0; color=colours[k], label=@sprintf("analytical, t = %.0f h", t / 3600))
    scatter!(
        p,
        x[1:4:end] .* 100,
        c[1:4:end] ./ c0;
        label=@sprintf("EKRem, t = %.0f h", t / 3600),
        marker=markers[k],
        markersize=sizes[k],
        markercolor=RGBA(1, 1, 1, 0),
        markerstrokecolor=colours[k],
    )
end
savefig(p, joinpath(output, "1D_trace_migration.png"))
println("\nFigure saved in ", joinpath(output, "1D_trace_migration.png"))
