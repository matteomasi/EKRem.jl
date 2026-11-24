#=
EKRem example: 1D setup, with binary electrolyte at steady state, under fixed applied voltage.

The setup consists of a 1 cm column of a water-saturated porous medium, with porosity of 0.4, 
and tortuosity factor of 0.5. The anolyte and catholyte are filled with NaCl at 10 and 100 mol/m³, 
respectively. The applied voltages are 0.5 V at the anode, and 0 V at the cathode. 
The model is run until the steady state is reached. The result is compared with the analytical 
at the steady state:

    c(x) = c_L + (c_R − c_L) x / L
    φ(x) = φ_L + (φ_R − φ_L) ln(c(x) / c_L) / ln(c_R / c_L)

At steady state the flux of each ion, Jᵢ = −Dᵢ (dcᵢ/dx + zᵢ cᵢ dψ/dx) with 
ψ = F φ / (R T), is constant along the column. Electroneutrality requires 
c₊ = c₋ = c, so that:

    J₊ / D₊ = −(dc/dx + c dψ/dx),    J₋ / D₋ = −(dc/dx − c dψ/dx).

The sum of the two equations shows that dc/dx is constant, so the concentration is linear.
Their difference shows that c dψ/dx is constant; with a linear c, its integration gives
the logarithmic potential, whose constant is fixed by the boundary potentials. The
diffusion coefficients enter only through the constants Jᵢ / Dᵢ, so the steady state
depends neither on them nor on the porosity and the tortuosity factor. The comparison is
repeated on three grids to show the convergence (second-order) of the discretisation.

Run from the repository root:
    julia --project=examples examples/1D_binary_fixed_potential.jl
=#

using EKRem
using Plots
using Printf

output = joinpath(@__DIR__, "output")
mkpath(output)

L = 0.01                   # column length (m)
cL, cR = 10.0, 100.0       # NaCl concentrations of the reservoirs (mol/m³)
φL, φR = 0.5, 0.0          # potentials of the reservoirs (V)
t_end = 1.0e6              # s; much longer than the diffusion time L² / (τ D) ≈ 1.5e5 s

electrolyte(c) = Dict("Na+" => c, "Cl-" => c)
species = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9)]

c_exact(x) = cL + (cR - cL) * x / L
φ_exact(x) = φL + (φR - φL) * log(c_exact(x) / cL) / log(cR / cL)

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

# Order of convergence p with respect to the grid spacing, assuming error ∝ hᵖ, estimated
# from two successive grids (the spacing is not exactly halved, hence the general formula)
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
p1 = plot(xf .* 1e3, c_exact.(xf); label="analytical", color=:black, xlabel="x (mm)", ylabel="c (mol/m³)")
scatter!(
    p1,
    x[1:5:end] .* 1e3,
    concentration(model, "Na+", t_end)[1:5:end];
    label="Na⁺ (EKRem)",
    markersize=5,
    marker=:square,
    markercolor=:grey,
)
scatter!(
    p1,
    x[1:5:end] .* 1e3,
    concentration(model, "Cl-", t_end)[1:5:end];
    label="Cl⁻ (EKRem)",
    markersize=4,
    marker=:circle,
    markercolor=:black,
)
p2 = plot(xf .* 1e3, φ_exact.(xf); label="analytical", color=:black, xlabel="x (mm)", ylabel="φ (V)")
scatter!(
    p2,
    x[1:5:end] .* 1e3,
    electric_potential(model, t_end)[1:5:end];
    label="EKRem",
    markersize=5,
    marker=:circle,
    markercolor=:white,
)
figure = plot(p1, p2; layout=(1, 2), size=(900, 350), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_binary_fixed_potential.png"))
println("Figure saved in ", joinpath(output, "1D_binary_fixed_potential.png"))
