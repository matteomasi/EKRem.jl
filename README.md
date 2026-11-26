# EKRem.jl

EKRem is a Julia package for the simulation and design of electrokinetic remediation systems. Multicomponent transport is described by the Nernst–Planck equations, the electric potential by the Poisson or the charge-conservation equation, and geochemical reactions are computed by [PhreeqcRM](https://github.com/usgs-coupled/phreeqcrm). The transport equations are discretised with finite volumes and Scharfetter–Gummel fluxes ([VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl)).

Three models are available, and they are organised by spatial dimension: `Model1D` for columns and laboratory setups, and `Model2D` and `Model3D` for field applications.

**Status:** pre-alpha. `Model1D` solves the Nernst–Planck equations with fixed concentrations at both ends. The electric potential is described by the Poisson equation (`Poisson()`), by charge conservation in an electroneutral solution (`ChargeConservation()`), or by a fixed field (`FixedField()`). The potential is fixed at both ends (voltage control), or at one end with a current density applied at the other (current control, with `ChargeConservation()`). Geochemical reactions are computed by PhreeqcRM at every node and coupled with transport by the sequential non-iterative approach (`database`, `input` and `coupling_step`; see the docstrings of `Model1D` and `run!`). Electrode compartments (`Compartment`) are well-mixed solutions with water electrolysis, under current or voltage control. Electroosmosis is described by a uniform electroosmotic permeability (`electroosmotic_permeability`). A zeta potential that depends on pH, `Model2D` and `Model3D` are not implemented yet.

## Example

Electromigration of a trace cation through a NaCl solution under an applied voltage of 1 V across a 10 cm column:

```julia
using EKRem

species = [Species("Na+", 1, 1.33e-9), Species("Cl-", -1, 2.03e-9), Species("T+", 1, 1.0e-9)]
background = Dict("Na+" => 100.0, "Cl-" => 100.0) # mol/m³

model = Model1D(;
    length=0.1, # m
    grid=(h_min=1.0e-5, h_max=1.0e-4, growth=1.1),
    species,
    initial=background,
    left=Reservoir(Dict("Na+" => 100.0, "Cl-" => 100.01, "T+" => 0.01); potential=1.0),
    right=Reservoir(background; potential=0.0),
    porosity=0.4,
    tortuosity=0.5,
)
run!(model, [1.0e4, 1.0e5]) # s

x = coordinates(model)
c = concentration(model, "T+", 1.0e5)
φ = electric_potential(model, 1.0e5)
i = current_density(model, 1.0e5) # A/m², between the nodes
```

For current control, `potential=1.0` of the left reservoir is replaced by `current_density=5.0` (A/m²), and `potential=ChargeConservation()` is passed to `Model1D`.

## Examples

Complete simulations are in the [`examples`](examples/) folder, in two groups; the file names begin with the dimension of the model, and each script saves a figure in the `output/` folder of its group.

Benchmarks (`examples/benchmarks/1D/`), compared with analytical solutions or with PHREEQC:

- `1D_trace_electromigration.jl`: demonstrates the coupled effects of electromigration and electroosmosis on the transport of a trace cation under a uniform electric field. The simulated concentration profiles are compared against the Ogata–Banks analytical solution (Ogata and Banks, 1961);
- `1D_binary_fixed_potential.jl`: demonstrates steady-state NaCl transport in a porous column under a fixed potential difference (0.5 V), comparing concentration and potential profiles with the analytical solution.
- `1D_binary_fixed_current.jl`: same as the previous example, under a fixed current density.
- `1D_benchmark_krabbenhoft.jl`: demonstrates migration of Na⁺, Cl⁻, K⁺, and OH⁻, solving the coupled Nernst–Planck–Poisson equations. Transient concentration profiles and steady-state distributions are compared with the analytical solution of Krabbenhøft and Krabbenhøft (2008);
- `1D_benchmark_phreeqc.jl`: validates the coupling between EKRem and PhreeqcRM. An acid front is generated at the anode and promotes nickel desorption from the porous medium. Results are compared against the Nernst–Planck module of PHREEQC (Appelo, 2017).

Examples (`examples/misc/`):

- `1D_compartments_fixed_voltage.jl`: simulates electrokinetic transport in a system with electrode compartments under a constant voltage difference, with aqueous speciation and water equilibrium computed by PhreeqcRM.

The examples have their own environment, which also contains the plotting package. They are run from the repository root, e.g.

```console
julia --project=examples examples/benchmarks/1D/1D_trace_electromigration.jl
```

or all together with `make examples`.

## Development set-up

The common development tasks are run with `make` (GNU Make); `make help` lists them:

```console
make install-dev    # instantiate the environments and install the development tools
make test           # run the unit tests (directly in the package environment)
make test-pkg       # run the unit tests with Pkg.test, in an isolated test environment
make examples       # run all examples
make docs           # build the documentation (site in docs/build/1)
make docs-rerun     # build the documentation, running all examples again
make docs-serve     # build the documentation and serve it locally with live reload
make lint           # check package quality with Aqua.jl
make format         # format the code with JuliaFormatter (settings in .JuliaFormatter.toml)
```

Without `make`, the environment is instantiated and the tests are run from the Julia package manager:

```julia
julia> ]
pkg> activate .
pkg> instantiate
pkg> test
```

## PhreeqcRM library

The PhreeqcRM shared library (version 3.8.6 or later) is not distributed with the package. Its path is set once, in the Julia REPL with the project activated:

```julia
using EKRem
set_phreeqcrm_library("C:/path/to/PhreeqcRM.dll")
```

The path is saved as a preference in `LocalPreferences.toml` next to the active `Project.toml`, and it is used automatically in later sessions and by the tests. Calling the function in the global environment makes the setting apply to all projects. `set_phreeqcrm_library(path; save=false)` loads a library for the current session only, and `phreeqcrm_library()` returns the path in use.

The library is selected in the following order of precedence: the library loaded in the current session, the environment variable `EKREM_PHREEQCRM`, and the saved preference.

Without a library, the PhreeqcRM tests are skipped. By default, the tests read the database from `../database/phreeqc.dat` relative to the library; another file can be given in `EKREM_PHREEQC_DATABASE`.


## References

- Appelo, C.A.J., 2017. Solute transport solved with the Nernst–Planck equation for concrete pores with ‘free’ water and a double layer. Cem. Concr. Res. 101, 102–113. https://doi.org/10.1016/j.cemconres.2017.08.030
- Krabbenhøft, K., Krabbenhøft, J., 2008. Application of the Poisson–Nernst–Planck equations to the migration test. Cem. Concr. Res. 38, 77–88. https://doi.org/10.1016/j.cemconres.2007.08.006
- Masi, M., Ceccarini, A., Iannelli, R., 2017. Multispecies reactive transport modelling of electrokinetic remediation of harbour sediments. J. Hazard. Mater. 326, 187–196. https://doi.org/10.1016/j.jhazmat.2016.12.032
- Masi, M., Paz-Garcia, J.M., Gomez-Lahoz, C., Villen-Guzman, M., Ceccarini, A., Iannelli, R., 2019. Modeling of electrokinetic remediation combining local chemical equilibrium and chemical reaction kinetics. J. Hazard. Mater. 371, 728–733. https://doi.org/10.1016/j.jhazmat.2019.03.014
- Ogata, A., Banks, R.B., 1961. A solution of the differential equation of longitudinal dispersion in porous media. U.S. Geological Survey Professional Paper 411-A, A1–A7. https://doi.org/10.3133/pp411A
- Rolle, M., Sprocati, R., Masi, M., Jin, B., Muniruzzaman, M., 2018. Nernst–Planck-based description of transport, Coulombic interactions, and geochemical reactions in porous media: Modeling approach and benchmark experiments. Water Resour. Res. 54, 3176–3195. https://doi.org/10.1002/2017WR022344
- Sprocati, R., Masi, M., Muniruzzaman, M., Rolle, M., 2019. Modeling electrokinetic transport and biogeochemical reactions in porous media: A multidimensional Nernst–Planck–Poisson approach with PHREEQC coupling. Adv. Water Resour. 127, 134–147. https://doi.org/10.1016/j.advwatres.2019.03.011

## License

Copyright (c) 2026 Matteo Masi

Licensed under the [Apache License, Version 2.0](LICENSE).
