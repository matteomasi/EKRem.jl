#=
# PHREEQC Ampere benchmark

The objective of this example is to showcase the coupling between SimEK and PhreeqcRM.
An acid front is generated at the anode, which desorbs nickel from the surface of a
porous medium. PhreeqcRM computes the surface complexation and the aqueous speciation. The
case is based on the benchmark of Sprocati et al. (2019), with different parameter values.
It is solved here in one dimension, and the results are compared with those of the 
Nernst–Planck module of PHREEQC version "Ampere" (Appelo, 2017), which can be found at:
https://www.hydrochemistry.eu/exmpls/electro_dif.html.

## Setup

The column has a length of 0.30 m and a cross-sectional area of 0.01 m². The porosity is 0.4,
and the tortuosity factor is 0.2. The pore water contains Na⁺ (10 mol/m³), Cl⁻ (12 mol/m³)
and nickel (1 mol/m³) at pH 7, and it is in equilibrium with the surface of the solid matrix.
A potential difference of 20 V is applied between the electrodes. The anode is placed in a closed
electrode compartment of 4 L. The cathode compartment is assumed to be completely buffered, 
so the concentrations at ``x = L`` are fixed.

Three surface reactions are defined (Sprocati et al., 2019, Eqs. 18–20):

```math
\mathrm{SurfOH + H^+ \rightleftharpoons SurfOH_2^+}, \qquad \log K = 3.8,
```
```math
\mathrm{SurfOH \rightleftharpoons SurfO^- + H^+}, \qquad \log K = -9.4,
```
```math
\mathrm{SurfOH + Ni^{2+} \rightleftharpoons SurfONi^+ + H^+}, \qquad \log K = -2.6.
```

The surface model has no electrostatic correction, and the surface has 0.015 mol of sites
per litre of pore water. The aqueous species and their equilibrium constants are taken from
the database `minteq.v4.dat`, distributed with PhreeqcRM. The diffusion coefficients of 
H⁺, OH⁻, Na⁺, Cl⁻, Ni²⁺ and NiCl⁺ are from Sprocati et al. (2019). The other aqueous species 
have a diffusion coefficient of ``10^{-9}`` m²/s.

## Reference solution

PHREEQC (Ampere) computes the reference solution with the input file
`data/1D_benchmark_phreeqc_ampere.phr`, on 240 cells with explicit time steps of 30 s. The
input file contains the same chemical definitions. The results are saved in the files
`data/1D_benchmark_phreeqc_profiles.tsv` and `data/1D_benchmark_phreeqc_current.tsv`. In
PHREEQC, the anode compartment is a boundary cell with 4 kg of water in equilibrium with
O₂(g), and the cathode is a boundary cell with a large volume of water. PHREEQC computes the
pore diffusion coefficient as ``D_w n^m``; the exponent ``m = \ln \tau / \ln n = 1.756``
gives the effective diffusion coefficient ``n \tau D_w`` of SimEK
(`-multi_d true 1e-9 0.4 0.0 1.756`).

The PhreeqcRM library must be set (see [Installation](../../../guide/installation.md)). The
example script can be run from the repository root:

```console
julia --project=examples examples/benchmarks/1D/1D_benchmark_phreeqc.jl
```
=#

using SimEK
using Plots
using Printf
default(; guidefontsize=13, tickfontsize=11, legendfontsize=10, linewidth=2)   # titles: default size

output = joinpath(@__DIR__, "output") #src
mkpath(output) #src

library = phreeqcrm_library()
library === nothing && error("The PhreeqcRM library is not set; see set_phreeqcrm_library")
database = joinpath(dirname(dirname(library)), "database", "minteq.v4.dat")
#md nothing #hide

# Reactions from minteq.v4 database, redefined with the diffusion coefficients (-dw). 
# The surface complexation are without electrostatic correction (-no_edl). 
# The amounts of the surface are per litre of pore water.
input = """
SOLUTION_SPECIES
H+ = H+;  log_k 0;  -dw 9.31e-9
Na+ = Na+;  log_k 0;  -dw 1.33e-9
Cl- = Cl-;  log_k 0;  -dw 2.03e-9
H2O = OH- + H+;  log_k -13.997;  delta_h 55.81 kJ;  -gamma 3.5 0;  -dw 5.27e-9
Ni+2 = Ni+2;  log_k 0;  -dw 6.6e-10
Ni+2 + Cl- = NiCl+;  log_k 0.408;  delta_h 2 kJ;  -gamma 0 0;  -dw 1e-9
SURFACE_MASTER_SPECIES
Surf SurfOH
SURFACE_SPECIES
SurfOH = SurfOH;  log_k 0
SurfOH + H+ = SurfOH2+;  log_k 3.8
SurfOH = SurfO- + H+;  log_k -9.4
SurfOH + Ni+2 = SurfONi+ + H+;  log_k -2.6
SOLUTION 1
pH 7;  Na 10;  Cl 12 charge;  Ni 1.0
SURFACE 1
SurfOH 0.015 24 1500
-equilibrate 1
-no_edl
SELECTED_OUTPUT 1
-reset false
-pH true
-molalities SurfONi+ SurfOH2+ SurfOH
END
"""
#md nothing #hide

# Definition of the transported species
# All aqueous species except water are selected. For the species without a diffusion coefficient
# in the database, a value of 1e-9 m²/s is used.
chemistry = Chemistry(1; database, input, initial=(solution=1,))
species = [Species(s.name, s.z, s.D > 0 ? s.D : 1.0e-9) for s in aqueous_species(chemistry) if s.name != "H2O"]
close(chemistry)

# The electric field reaches about 80 V/m, and the potential drop along an edge of length h
# is ``|z \Delta\psi| \approx 3 \times 10^3 h`` (h in m). Where this cell Péclet number is
# much larger than 1, the Scharfetter–Gummel flux approaches upwinding and smooths the
# concentration profiles. The grid is therefore refined toward the electrodes (0.25 mm).
L = 0.3                    # column length (m)
A = 0.01                   # cross-sectional area (m²)
porosity, tortuosity = 0.4, 0.2
model = Model1D(;
    length=L,
    grid=(h_min=2.5e-4, h_max=1.25e-3, growth=1.1),
    species,
    initial=(solution=1, surface=1),
    left=Compartment(; initial=(solution=1,), volume=4.0e-3, potential=20.0),   # anode
    right=Reservoir(; solution=1, potential=0.0),                               # cathode
    potential=ChargeConservation(),
    porosity,
    tortuosity,
    area=A,
    database,
    input,
)
display(model) #src
println() #src
#md model

# The model is run for 80 h, with results every 2 h and a coupling step of 300 s.
hour = 3600.0
times = collect(2.0:2.0:80.0) .* hour
t_end = times[end]
elapsed = @elapsed run!(model, times; coupling_step=300.0) #src
@printf("Run time: %.0f s\n", elapsed) #src
#md run!(model, times; coupling_step=300.0, progress=false)
#md nothing #hide

# The reference solution of PHREEQC (Ampere) contains the profiles at 80 h and the current
# every 2 h.
#hide-begin
# Each line of the files written by PHREEQC begins with a tab.
function read_table(file)
    lines = filter(!isempty, strip.(readlines(file)))
    header = split(lines[1], '\t')
    rows = [parse.(Float64, split(line, '\t')) for line in lines[2:end]]
    return Dict(String(h) => [row[k] for row in rows] for (k, h) in enumerate(header))
end
data = joinpath(pkgdir(SimEK), "examples", "benchmarks", "1D", "data")
ref = read_table(joinpath(data, "1D_benchmark_phreeqc_profiles.tsv"))
ref_current = read_table(joinpath(data, "1D_benchmark_phreeqc_current.tsv"))

J = [current_density(model, t)[1] for t in times]             # A/m²
J_ref = ref_current["current_A"] ./ A
@printf("Current density at 80 h: SimEK %.3f A/m², PHREEQC %.3f A/m²\n", J[end], J_ref[end]) #src
#hide-end
#md nothing #hide

# In the figure below SimEK (lines) is compared with PHREEQC (points).
# In the panel (a) the concentrations of the ions and pH are shown, together with the pH on the right axis.
# The panel (b) shows the surface species concentrations. The panel (c) shows the contribution 
# of each ion to the electrical conductivity of the porous medium, calculated as
# ``\sigma_i = F^2 z_i^2 n \tau D_i c_i / (R T)``.
# The panel (d) shows the evolution of the current density over time.
#hide-begin
σ(z, D, c) = FARADAY^2 / (GAS_CONSTANT * 298.15) * z^2 * porosity * tortuosity * D .* c

# Profiles at 80 h: SimEK (lines) and PHREEQC (points, in every sixth cell)
x = coordinates(model)
x_ref = ref["x_m"]
marks = 1:6:length(x_ref)
so = selected_output(model, t_end)
ions = [   # species, column of the reference, label, scale factor, colour
    ("Ni+2", "Ni", "Ni²⁺", 1, "#0072B2"),
    ("NiCl+", "NiCl", "NiCl⁺ × 10", 10, "#CC79A7"),
    ("Na+", "Na", "Na⁺", 1, "#009E73"),
    ("Cl-", "Cl", "Cl⁻", 1, "#D55E00"),
]
D = Dict(s.name => s.D for s in species)
z = Dict(s.name => s.z for s in species)
marker = (; markerstrokewidth=0, markersize=3.5)

# (a) Concentrations and pH (right axis); the legend lists both axes and lies in the empty
# region between the profiles
p1 = plot(; title="(a) Concentration", xlabel="x (cm)", ylabel="c (mol/m³)", ylims=(0, 36), legend=(0.55, 0.72))
for (name, col, label, f, color) in ions
    plot!(p1, x .* 100, concentration(model, name, t_end) .* f; color, label)
    scatter!(p1, x_ref[marks] .* 100, ref[col][marks] .* f; color, label="", marker...)
end
plot!(p1, [NaN], [NaN]; color=:black, label="pH (right axis)")
scatter!(p1, [NaN], [NaN]; color=:gray, label="PHREEQC", marker...)
p1_pH = twinx(p1)
plot!(p1_pH, x .* 100, so["pH"]; color=:black, ylims=(0, 8), ylabel="pH", legend=false)
scatter!(p1_pH, x_ref[marks] .* 100, ref["pH"][marks]; color=:black, marker...)

# (b) Surface species (without the boundary nodes and cells, which have no surface)
p2 = plot(; title="(b) Surface species", xlabel="x (cm)", ylabel="c (mol/kgw)", legend=:right)
inner = 2:(length(x) - 1)
inner_marks = filter(k -> 1 < k < length(x_ref), marks)
surface_species = [
    ("SurfONi+", "SurfONi", "SurfONi⁺", "#0072B2"),
    ("SurfOH2+", "SurfOH2", "SurfOH₂⁺", "#D55E00"),
    ("SurfOH", "SurfOH", "SurfOH", "#009E73"),
]
for (name, col, label, color) in surface_species
    plot!(p2, x[inner] .* 100, so["m_$name(mol/kgw)"][inner]; color, label)
    scatter!(p2, x_ref[inner_marks] .* 100, ref[col][inner_marks]; color, label="", marker...)
end
scatter!(p2, [NaN], [NaN]; color=:gray, label="PHREEQC", marker...)

# (c) Contribution of each ion to the conductivity
p3 = plot(; title="(c) Conductivity", xlabel="x (cm)", ylabel="σᵢ (S/m)", legend=:topright)
for (name, col, label, f, color) in [ions; ("H+", "H", "H⁺", 1, "#E69F00")]
    plot!(p3, x .* 100, σ(z[name], D[name], concentration(model, name, t_end)) .* f; color, label)
    scatter!(p3, x_ref[marks] .* 100, σ(z[name], D[name], ref[col][marks]) .* f; color, label="", marker...)
end
scatter!(p3, [NaN], [NaN]; color=:gray, label="PHREEQC", marker...)

# (d) Current density as a function of time
p4 = plot(; title="(d) Current density", xlabel="t (h)", ylabel="J (A/m²)", legend=:bottomright)
plot!(p4, times ./ hour, J; color=:black, label="SimEK")
scatter!(p4, ref_current["time_h"][2:end], J_ref[2:end]; color=:black, label="PHREEQC", marker...)

figure = plot(p1, p2, p3, p4; layout=(2, 2), size=(1000, 800), margin=5Plots.mm)
savefig(figure, joinpath(output, "1D_benchmark_phreeqc.png")) #src
println("Figure saved in ", joinpath(output, "1D_benchmark_phreeqc.png")) #src
#hide-end
#md figure #hide

#=
## References
Appelo, C.A.J., 2017. Solute transport solved with the Nernst–Planck equation for concrete 
pores with ‘free’ water and a double layer. Cem. Concr. Res. 101, 102–113. 
https://doi.org/10.1016/J.CEMCONRES.2017.08.030 .

Sprocati, R., Masi, M., Muniruzzaman, M., Rolle, M., 2019. Modeling electrokinetic transport
and biogeochemical reactions in porous media: A multidimensional Nernst–Planck–Poisson
approach with PHREEQC coupling. Adv. Water Resour. 127, 134–147.
https://doi.org/10.1016/j.advwatres.2019.03.011
=#
