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

Complete simulations are in the [`examples`](examples/) folder; the file names begin with the dimension of the model. Each example compares the results with an analytical solution and saves a figure in `examples/output/`:

- `1D_trace_electromigration.jl`: electromigration and electroosmosis of a trace ion in a uniform electric field (Ogata–Banks solution);
- `1D_binary_fixed_potential.jl`: steady state of a binary electrolyte between two reservoirs at different concentrations and potentials, on three grids;
- `1D_binary_fixed_current.jl`: the same column under a fixed current density (charge-conservation formulation): potential at the anode as a function of time, and steady state;
- `1D_benchmark_krabbenhoft.jl`: migration test with four monovalent ions (Na⁺, Cl⁻, K⁺, OH⁻) under a fixed electric potential difference, Nernst–Planck–Poisson equations: transient profiles, and steady state compared with the analytical solution of Krabbenhøft and Krabbenhøft (2008).

The examples have their own environment, which also contains the plotting package. They are run from the repository root, e.g.

```console
julia --project=examples examples/1D_trace_electromigration.jl
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
