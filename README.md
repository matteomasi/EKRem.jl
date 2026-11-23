# EKRem.jl

EKRem is a Julia package for the simulation and design of electrokinetic remediation systems. Multicomponent transport is described by the Nernst–Planck equations, the electric potential by the Poisson or the charge-conservation equation, and geochemical reactions are computed by [PhreeqcRM](https://github.com/usgs-coupled/phreeqcrm). The transport equations are discretised with finite volumes and Scharfetter–Gummel fluxes ([VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl)).

Three models are available, and they are organised by spatial dimension: `Model1D` for columns and laboratory setups, and `Model2D` and `Model3D` for field applications.

**Status:** pre-alpha. Only the package skeleton, the Scharfetter–Gummel flux and the first functions of the PhreeqcRM interface are available.

## Development set-up

The common development tasks are run with `make` (GNU Make); `make help` lists them:

```console
make install-dev    # instantiate the environment and install the development tools
make test           # run the test suite
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

The PhreeqcRM shared library is not distributed with the package. Its path is set once, in the Julia REPL with the project activated:

```julia
using EKRem
set_phreeqcrm_library("C:/path/to/PhreeqcRM.dll")
```

The path is saved as a preference in `LocalPreferences.toml` next to the active `Project.toml`, and it is used automatically in later sessions and by the tests. Calling the function in the global environment makes the setting apply to all projects. `set_phreeqcrm_library(path; save=false)` loads a library for the current session only, and `phreeqcrm_library()` returns the path in use.

The library is selected in the following order of precedence: the library loaded in the current session, the environment variable `EKREM_PHREEQCRM`, and the saved preference.

Without a library, the PhreeqcRM tests are skipped. By default, the tests read the database from `../database/phreeqc.dat` relative to the library; another file can be given in `EKREM_PHREEQC_DATABASE`.
