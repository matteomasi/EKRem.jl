#=
EKRem example: 1D setup, with binary electrolyte, under a fixed applied current

The setup consists of a 1 cm column of a water-saturated porous medium, with porosity of 0.4, 
and tortuosity factor of 0.5. The anolyte and catholyte are filled with NaCl at 10 and 100 mol/m³, 
respectively.
A fixed current density of 5 A/m² is applied at the anode, while the potential at the cathode is 
fixed at 0 V. The charge-conservation equation is used to describe the potential. 
The simulation is computed until the steady state is reached, and compared with
the analytical solution.

At steady state the flux of each ion, Jᵢ = −Dᵢ* (dcᵢ/dx + zᵢ cᵢ dψ/dx) with ψ = F φ / (R T)
and Dᵢ* = n τ Dᵢ, is constant along the column. Electroneutrality requires c₊ = c₋ = c, so
that

    J₊ / D₊* + J₋ / D₋* = −2 dc/dx,    J₊ / D₊* − J₋ / D₋* = −2 c dψ/dx.

The first equation shows that the concentration is linear, c(x) = c_L + (c_R − c_L) x / L.
The current density i = F (J₊ − J₋) is fixed, so that the first equation and
J₊ − J₋ = i / F determine both fluxes:

    J₋ = −(2 dc/dx + i / (F D₊*)) / (1/D₊* + 1/D₋*),    J₊ = J₋ + i / F.

The second equation then gives c dψ/dx = B, with B = −(J₊ / D₊* − J₋ / D₋*) / 2, and its
integration with φ(L) = 0 gives

    φ(x) = (R T / F) (B / (dc/dx)) ln(c(x) / c_R).

Unlike the steady state under a fixed voltage, this solution depends on the diffusion
coefficients, the porosity and the tortuosity factor, through the conductivity of the
column.

Run from the repository root:
    julia --project=examples examples/1D_binary_fixed_current.jl
=#

using EKRem
using Plots
using Printf

output = joinpath(@__DIR__, "output")
mkpath(output)

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
display(model)
println()

# Output times from 10 s to 10⁶ s, much longer than the diffusion time L² / (τ D) ≈ 1.5e5 s
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
p1 = plot(xf .* 1e3, c_exact.(xf); label="analytical", color=:black, xlabel="x (mm)", ylabel="c (mol/m³)")
scatter!(
    p1,
    x[1:10:end] .* 1e3,
    concentration(model, "Na+", t_end)[1:10:end];
    label="EKRem (steady state)",
    markercolor=:white,
    markersize=5,
)
p2 = plot(xf .* 1e3, φ_exact.(xf); label="analytical", color=:black, xlabel="x (mm)", ylabel="φ (V)")
scatter!(
    p2,
    x[1:10:end] .* 1e3,
    electric_potential(model, t_end)[1:10:end];
    label="EKRem (steady state)",
    markercolor=:white,
    markersize=5,
)
p3 = plot(
    output_times(model)[2:end],
    φ_an[2:end];
    xscale=:log10,
    color=:black,
    label="EKRem",
    xlabel="t (s)",
    ylabel="φ anode (V)",
)
hline!(p3, [φ_exact(0.0)]; color=:black, linestyle=:dash, label="steady state (analytical)")
figure = plot(p1, p2, p3; layout=(1, 3), size=(1200, 350), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_binary_fixed_current.png"))
println("Figure saved in ", joinpath(output, "1D_binary_fixed_current.png"))
