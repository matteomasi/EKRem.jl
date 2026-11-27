# SimEK.jl

SimEK is a Julia package for the simulation and design of electrokinetic (EK) remediation systems. 

Electrokinetic remediation (EKR) is a technique for removing contaminants from porous media by applying an electric field, which induces the transport of contaminants through the medium (Reddy and Cameselle, 2009). EKR has been successfully employed for the treatment of soils, sediments and groundwater contaminated by a wide variety of pollutants such as heavy metals, organic compounds, and radionuclides (Masi et al., 2017; Wen et al., 2021).

The mathematical description of EKR therefore requires the coupled simulation of species transport, electrochemical processes, and geochemical reactions (Sprocati et al., 2019). In SimEK, multicomponent electrokinetic transport is described by the Nernst–Planck equations, the electric potential by either the Poisson or the charge-conservation equation, and geochemical reactions are computed by [PhreeqcRM](https://github.com/usgs-coupled/phreeqcrm). 

The equations are discretised in space by the Voronoi finite volume method ([VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl)). The Nernst–Planck fluxes are approximated using the Scharfetter–Gummel scheme. The resulting differential–algebraic system is integrated in time with adaptive implicit methods from [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl).

SimEK is intended to support the interpretation of laboratory experiments and the design of field treatments.
Depending on the geometry of the system and the dimensionality of the electric field, EKR configurations can be represented using one-, two-, or three-dimensional models. Accordingly, SimEK provides three model classes:
- `Model1D` for columns and laboratory cells;
- `Model2D` for two-dimensional laboratory setups and cross-sections of field installations;
- `Model3D` for full three-dimensional field installations.


## Status

SimEK is in a pre-alpha stage; currently, only `Model1D` is implemented. It includes:

- the electric potential described by the Poisson equation (`Poisson()`), by charge conservation in an electroneutral solution (`ChargeConservation()`) or by a fixed field (`FixedField()`);
- boundaries in contact with reservoirs of fixed composition (`Reservoir`) or with well-mixed electrode compartments in which water electrolysis takes place (`Compartment`);
- controlled voltage, or controlled current;
- electroosmosis with uniform electroosmotic permeability;
- geochemical reactions computed by PhreeqcRM at every node and coupled with transport by the sequential non-iterative approach.

`Model2D`, `Model3D` are not implemented yet.

## Getting started

SimEK requires Julia 1.10 or later. The package is not registered yet; it is installed from this repository with the Julia package manager:

```julia
julia> ]
pkg> add https://github.com/matteomasi/SimEK.jl
```

Models with geochemical reactions also require the PhreeqcRM shared library, which is installed separately (see [PhreeqcRM library](#phreeqcrm-library)). Models without geochemical reactions, such as the example below, do not require it.

## Basic usage 

The following example simulates the electromigration of a trace cation (T⁺) through a 10 cm porous column saturated with a 0.1 M NaCl solution. The cation enters from the reservoir at the anode, and a voltage of 1 V is applied between the two ends of the column. The potential is described by the Poisson equation, and the profiles are saved after 10⁴ and 10⁵ s.

```julia
using SimEK

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

x = coordinates(model) # Get the coordinates of the grid nodes
c = concentration(model, "T+", 1.0e5) # Get the concentration of T⁺ at 10⁵ s
φ = electric_potential(model, 1.0e5) # Get the electric potential at 10⁵ s
J = current_density(model, 1.0e5) # Get the current density (A/m²) at 10⁵ s
```

## Examples

Complete simulations are available in the [`examples`](examples/) folder, including benchmarks against analytical solutions and PHREEQC. They are described in the [examples section](docs/src/examples.md) of the documentation.

## Development

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
using SimEK
set_phreeqcrm_library("C:/path/to/PhreeqcRM.dll")
```

The path is saved as a preference in `LocalPreferences.toml` next to the active `Project.toml`, and it is used automatically in later sessions and by the tests. Calling the function in the global environment makes the setting apply to all projects. `set_phreeqcrm_library(path; save=false)` loads a library for the current session only, and `phreeqcrm_library()` returns the path in use.

The library is selected in the following order of precedence: the library loaded in the current session, the environment variable `SIMEK_PHREEQCRM`, and the saved preference.

Without a library, the PhreeqcRM tests are skipped. By default, the tests read the database from `../database/phreeqc.dat` relative to the library; another file can be given in `SIMEK_PHREEQC_DATABASE`.


## Contributing

Contributions to SimEK are welcome, from questions and bug reports to pull requests with new features, fixes or documentation. Questions and bug reports are submitted through the [issue tracker](https://github.com/matteomasi/SimEK.jl/issues). Changes to the code should be discussed in an issue before a pull request is opened. The procedure and the conventions for tests, examples and documentation are described in [CONTRIBUTING.md](CONTRIBUTING.md).


## License

Copyright (c) 2026 Matteo Masi

Licensed under the [Apache License, Version 2.0](LICENSE).


## References

- Reddy, K.R., Cameselle, C., 2009. Electrochemical Remediation Technologies for Polluted Soils, Sediments and Groundwater. John Wiley & Sons, Hoboken, NJ. https://doi.org/10.1002/9780470523650

- Masi, M., Ceccarini, A., Iannelli, R., 2017. Multispecies reactive transport modelling of electrokinetic remediation of harbour sediments. J. Hazard. Mater. 326, 187–196. https://doi.org/10.1016/j.jhazmat.2016.12.032

- Wen, D., Fu, R., Li, Q., 2021. Removal of inorganic contaminants in soil by electrokinetic remediation technologies: A review. J. Hazard. Mater. 401, 123345. https://doi.org/10.1016/j.jhazmat.2020.123345

- Sprocati, R., Masi, M., Muniruzzaman, M., Rolle, M., 2019. Modeling electrokinetic transport and biogeochemical reactions in porous media: A multidimensional Nernst–Planck–Poisson approach with PHREEQC coupling. Adv. Water Resour. 127, 134–147. https://doi.org/10.1016/j.advwatres.2019.03.011

