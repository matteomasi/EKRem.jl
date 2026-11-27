# Changelog

All notable changes to SimEK are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). Before
version 1.0.0, a minor version may contain breaking changes.

## [Unreleased]

### Changed

- The package has been renamed from EKRem to SimEK. The environment variables are now
  `SIMEK_PHREEQCRM` and `SIMEK_PHREEQC_DATABASE`, and the saved preference is read from
  the `[SimEK]` section of `LocalPreferences.toml`.

### Added

- `Model1D`: one-dimensional multicomponent transport described by the Nernst–Planck
  equations, discretised with finite volumes and Scharfetter–Gummel fluxes (VoronoiFVM.jl),
  on uniform or refined grids.
- Description of the electric potential by the Poisson equation (`Poisson()`), by charge
  conservation in an electroneutral solution (`ChargeConservation()`) or by a fixed field
  (`FixedField()`).
- Voltage control (potential fixed at both ends) and current control (potential fixed at
  one end and current density applied at the other, with `ChargeConservation()`).
- Reservoir boundary conditions with fixed concentrations (`Reservoir`).
- Electrode compartments (`Compartment`): well-mixed solutions with water electrolysis,
  under current or voltage control.
- Electroosmosis described by a uniform electroosmotic permeability.
- Interface to the C API of PhreeqcRM, with the path of the shared library set by
  `set_phreeqcrm_library`, by the environment variable `SIMEK_PHREEQCRM` or by a saved
  preference.
- Coupling of transport and geochemical reactions by the sequential non-iterative approach.
- Benchmarks against analytical solutions (Ogata–Banks; Krabbenhøft and Krabbenhøft, 2008)
  and against the Nernst–Planck module of PHREEQC, and an example with electrode
  compartments.
- Documentation built with Documenter.jl and DocumenterVitepress.jl, with the examples
  rendered by Literate.jl.
- Continuous integration of the tests on Linux, Windows and macOS, with PhreeqcRM built from
  its source distribution.

[Unreleased]: https://github.com/matteomasi/SimEK.jl/commits/main
