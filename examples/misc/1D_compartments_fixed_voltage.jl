#=
# Electrokinetic transport with electrode compartments under controlled voltage

The setup consists of a 10 cm water-saturated porous column (porosity 0.4, tortuosity factor 0.5, 
cross-sectional area 50 cm²), initially containing 10 mol/m³ NaCl, and two 0.5 L electrode compartments 
filled with the same solution, one at each end of the column. 
A voltage of 10 V is applied at the anode, driving transport and water electrolysis.
Under voltage-controlled conditions, the electrolysis rate is governed by the electric current 
through the column, which is calculated by the model (see [`Compartment`](@ref)). 
The resulting acidic and alkaline fronts migrate into the porous medium and neutralise upon mixing. 
Aqueous speciation, including water dissociation, is computed using PhreeqcRM.

The PhreeqcRM library must be set (see [Installation](../../guide/installation.md)); the
database `phreeqc.dat` is taken from the `database` folder of the PhreeqcRM installation.

The example script can be run from the repository root:

```console
julia --project=examples examples/misc/1D_compartments_fixed_voltage.jl
```
=#

using SimEK
using Plots
using Printf
default(; guidefontsize=13, tickfontsize=11, legendfontsize=10)   # titles: default size

output = joinpath(@__DIR__, "output") #src
mkpath(output) #src

library = phreeqcrm_library()
library === nothing && error("The PhreeqcRM library is not set; see set_phreeqcrm_library")
database = joinpath(dirname(dirname(library)), "database", "phreeqc.dat")

input = """
SOLUTION 1
    pH 7 charge
    Na 10
    Cl 10
SELECTED_OUTPUT 1
    -reset false
    -pH true
END
"""

L = 0.1                    # column length (m)
U = 10.0                   # applied voltage (V)
compartment(potential) = Compartment(; initial=(solution=1,), volume=5.0e-4, potential)
model = Model1D(;
    length=L,
    grid=collect(range(0, L; length=51)),
    initial=(solution=1,),
    left=compartment(U),         # anode
    right=compartment(0.0),      # cathode
    potential=ChargeConservation(),
    porosity=0.4,
    tortuosity=0.5,
    area=5.0e-3,
    database,
    input,
)
display(model) #src
println() #src
#md model

# The model is run for 2 days, with results every hour
hour = 3600.0
times = collect(1.0:48.0) .* hour
run!(model, times; coupling_step=300.0)

# Current density J and pH of the compartments (first and last node) as a function of time
t = output_times(model)
J = [current_density(model, tk)[1] for tk in t]
pH_anode = [selected_output(model, tk)["pH"][1] for tk in t]
pH_cathode = [selected_output(model, tk)["pH"][end] for tk in t]
charge = sum((J[1:(end - 1)] .+ J[2:end]) ./ 2 .* diff(t)) / FARADAY   # trapezoidal rule
@printf(
    "Current density J: %.2f A/m² at t = 0, %.2f A/m² at t = 48 h
",
    J[1],
    J[end]
)
@printf(
    "Charge passed in 48 h: %.3f mol/m² (H⁺ produced at the anode, OH⁻ at the cathode)
",
    charge
)
@printf(
    "pH of the compartments at 48 h: anode %.2f, cathode %.2f
",
    pH_anode[end],
    pH_cathode[end]
)

# Profiles in the column (the compartment nodes are excluded)
x = coordinates(model)[2:(end - 1)] .* 100
column(values) = values[2:(end - 1)]
profile_times = (1.0, 2.0, 6.0, 48.0) .* hour
colors = ["#0072B2", "#D55E00", "#009E73", "#000000"]   # Okabe–Ito palette (colour-blind safe)
p1 = plot(; title="pH", xlabel="x (cm)", ylabel="pH", legend=:topleft)
p2 = plot(; title="Na⁺ (solid) and Cl⁻ (dashed)", xlabel="x (cm)", ylabel="c (mol/m³)", legend=false)
for (k, tk) in enumerate(profile_times)
    label = @sprintf("t = %.0f h", tk / hour)
    plot!(p1, x, column(selected_output(model, tk)["pH"]); color=colors[k], linewidth=2, label)
    plot!(p2, x, column(concentration(model, "Na+", tk)); color=colors[k], linewidth=2)
    plot!(p2, x, column(concentration(model, "Cl-", tk)); color=colors[k], linewidth=2, linestyle=:dash)
end
p3 = plot(; title="Electrolytes", xlabel="t (h)", ylabel="pH", legend=:right)
plot!(p3, t ./ hour, pH_anode; color=colors[2], linewidth=2, label="anode chamber")
plot!(p3, t ./ hour, pH_cathode; color=colors[1], linewidth=2, label="cathode chamber")
p4 = plot(; title="Current density", xlabel="t (h)", ylabel="J (A/m²)", legend=false)
plot!(p4, t ./ hour, J; color=:black, linewidth=2)
figure = plot(p1, p2, p3, p4; layout=(2, 2), size=(1000, 800), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_compartments_fixed_voltage.png")) #src
println("Figure saved in ", joinpath(output, "1D_compartments_fixed_voltage.png")) #src
#md figure
