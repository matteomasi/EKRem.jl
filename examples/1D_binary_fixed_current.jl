#=
# Binary electrolyte under a fixed current

A 1 cm column of a water-saturated porous medium (porosity 0.4, tortuosity factor 0.5) lies
between two reservoirs with NaCl at 10 mol/m³ (anolyte) and 100 mol/m³ (catholyte). A
current density of 5 A/m² is applied at the anode, and the potential at the cathode is
fixed at 0 V. The potential is described by the charge-conservation equation. The model is
run until the steady state is reached, and the result is compared with the analytical
solution.

At steady state, the flux of each ion,
``J_i = -D_i^* \left( \mathrm{d}c_i/\mathrm{d}x + z_i c_i \, \mathrm{d}\psi/\mathrm{d}x \right)``
with ``\psi = F \phi / (R T)`` and ``D_i^* = n \tau D_i``, is constant along the column.
Electroneutrality requires ``c_+ = c_- = c``, so that

```math
\frac{J_+}{D_+^*} + \frac{J_-}{D_-^*} = -2 \frac{\mathrm{d}c}{\mathrm{d}x},
\qquad
\frac{J_+}{D_+^*} - \frac{J_-}{D_-^*} = -2 c \frac{\mathrm{d}\psi}{\mathrm{d}x}.
```

The first equation shows that the concentration is linear,
``c(x) = c_L + (c_R - c_L) x / L``. The current density ``i = F (J_+ - J_-)`` is fixed, so
that the first equation and ``J_+ - J_- = i / F`` determine both fluxes:

```math
J_- = -\frac{2 \, \mathrm{d}c/\mathrm{d}x + i / (F D_+^*)}{1/D_+^* + 1/D_-^*},
\qquad J_+ = J_- + \frac{i}{F}.
```

The second equation then gives ``c \, \mathrm{d}\psi/\mathrm{d}x = B``, with
``B = -(J_+ / D_+^* - J_- / D_-^*) / 2``, and its integration with ``\phi(L) = 0`` gives

```math
\phi(x) = \frac{R T}{F} \frac{B}{\mathrm{d}c/\mathrm{d}x} \ln\frac{c(x)}{c_R}.
```

Unlike the steady state under a fixed voltage, this solution depends on the diffusion
coefficients, the porosity and the tortuosity factor, through the conductivity of the
column.

Run from the repository root:

```console
julia --project=examples examples/1D_binary_fixed_current.jl
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
i = 5.0                    # current density applied at the left boundary (A/m²)
n, τ = 0.4, 0.5            # porosity and tortuosity factor
T = 298.15                 # temperature (K)
D₊, D₋ = 1.33e-9, 2.03e-9  # diffusion coefficients of Na⁺ and Cl⁻ in water (m²/s)

electrolyte(c) = Dict("Na+" => c, "Cl-" => c)
species = [Species("Na+", 1, D₊), Species("Cl-", -1, D₋)]

# Analytical steady state
Dp, Dm = n * τ * D₊, n * τ * D₋
dcdx = (cR - cL) / L
Jm = -(2 * dcdx + i / (FARADAY * Dp)) / (1 / Dp + 1 / Dm)
Jp = Jm + i / FARADAY
B = -(Jp / Dp - Jm / Dm) / 2
c_exact(x) = cL + dcdx * x
φ_exact(x) = GAS_CONSTANT * T / FARADAY * B / dcdx * log(c_exact(x) / cR)

model = Model1D(;
    length=L,
    grid=collect(range(0, L; length=201)),
    species,
    initial=electrolyte(cL),
    left=Reservoir(electrolyte(cL); current_density=i),
    right=Reservoir(electrolyte(cR); potential=0.0),
    potential=ChargeConservation(),
    porosity=n,
    tortuosity=τ,
    temperature=T,
)
display(model) #src
println() #src
#md model

# Output times from 10 s to 10⁶ s, much longer than the diffusion time
# ``L^2 / (\tau D) \approx 1.5 \times 10^5`` s
times = 10.0 .^ range(1, 6; length=26)
println("\nSimulation")
run!(model, times)
t_end = times[end]

x = coordinates(model)
φ_an = [electric_potential(model, t)[1] for t in output_times(model)]   # potential at the anode
error_c = maximum(abs.(concentration(model, "Na+", t_end) .- c_exact.(x))) / cR
error_φ = maximum(abs.(electric_potential(model, t_end) .- φ_exact.(x))) / φ_exact(0.0)
i_model = current_density(model, t_end)

println("\nResults at t = $(t_end) s")
@printf("  potential at the anode φ_an: %.5f V (analytical steady state %.5f V)\n", φ_an[end], φ_exact(0.0))
@printf("  max error of c / c_R: %.2e;  max error of φ / φ_an: %.2e\n", error_c, error_φ)
@printf("  current density: %.6f to %.6f A/m² (applied %.1f A/m²)\n", extrema(i_model)..., i)

# Steady-state profiles and potential at the anode as a function of time
xf = range(0, L; length=400)
marks = 1:10:length(x)
p1 = plot(; title="Concentration (steady state)", xlabel="x (mm)", ylabel="c (mol/m³)", legend=:topleft)
scatter!(
    p1, x[marks] .* 1e3, concentration(model, "Na+", t_end)[marks]; label="EKRem", markercolor=:black, markersize=6
)
plot!(p1, xf .* 1e3, c_exact.(xf); color=:black, label="analytical")
p2 = plot(; title="Electric potential (steady state)", xlabel="x (mm)", ylabel="φ (V)", legend=:topright)
scatter!(p2, x[marks] .* 1e3, electric_potential(model, t_end)[marks]; label="EKRem", markercolor=:black, markersize=6)
plot!(p2, xf .* 1e3, φ_exact.(xf); color=:black, label="analytical")
p3 = plot(; title="Potential at the anode", xlabel="t (s)", ylabel="φ (V)", xscale=:log10, legend=:topright)
scatter!(p3, output_times(model)[2:end], φ_an[2:end]; label="EKRem", markercolor=:black, markersize=6)
hline!(p3, [φ_exact(0.0)]; color=:black, label="analytical (steady state)")
figure = plot(p1, p2, p3; layout=@layout([a b; c _]), size=(1000, 800), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_binary_fixed_current.png")) #src
println("Figure saved in ", joinpath(output, "1D_binary_fixed_current.png")) #src
#md figure
