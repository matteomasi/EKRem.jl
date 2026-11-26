# Contributing to EKRem

Contributions to EKRem are welcome, from questions and bug reports to pull requests with
new features, fixes or documentation. This document describes how each kind of
contribution is handled and how the project is maintained.

## Questions

1. Search the [issue tracker](https://github.com/matteomasi/EKRem/issues) for an existing
   issue on the same topic.
2. If none is found, open a new issue and apply the label `question`.

Questions on the theory (Nernst–Planck transport, electric potential, coupling with
PhreeqcRM) are also welcome, but the documentation (in `docs/`, built with `make docs`)
should be consulted first.

## Bug reports

1. Search the [issue tracker](https://github.com/matteomasi/EKRem/issues) for an existing
   issue on the same problem.
2. If none is found, open a new issue with the label `bug`. The report should contain:
   - a minimal script that reproduces the problem, and the full error message or a
     description of the unexpected result;
   - the version of EKRem (or the commit hash), the output of `versioninfo()` and the
     operating system;
   - for problems involving geochemical reactions, the version of PhreeqcRM, the
     thermodynamic database and the PHREEQC input.

Results that disagree with an analytical solution or with another code are reported as
bugs, together with the reference solution.

## Changes to the code

1. Open an issue that describes the proposed change *before starting the work*, and wait
   for a discussion of the approach. This avoids duplicated effort and changes that do not
   fit the design of the package.
2. Fork the repository and create a branch from the latest commit of `main`.
3. Set up the development environment (see the
   [development set-up](docs/src/development.md)):

   ```console
   make install-dev
   ```

   The tests of the geochemical coupling require the PhreeqcRM shared library, set with
   `set_phreeqcrm_library(path)` or with the environment variable `EKREM_PHREEQCRM`.
   Without the library, these tests are skipped.
4. Make the change, and add tests for it.
5. Check the change:

   ```console
   make test           # unit tests
   make lint           # package quality (Aqua.jl)
   make format         # formatting (JuliaFormatter, settings in .JuliaFormatter.toml)
   make docs           # documentation, if docstrings or pages were changed
   ```

6. Add an entry to the section `[Unreleased]` of [CHANGELOG.md](CHANGELOG.md).
7. Push the branch to the fork and open a pull request against `main`. The continuous
   integration runs the tests on Linux, Windows and macOS.

A contribution is not refused because tests or documentation are missing: the pull request
can be opened anyway, and help will be provided to complete it.

### Conventions

- **Tests.** The tests in `test/` are fast unit tests. They check results that follow from
  physics and chemistry (conservation of mass and charge, equilibrium, stoichiometry,
  analytical solutions), not values that may change between versions of a thermodynamic
  database or of a dependency. Complete simulations belong in `examples/`.
- **Examples.** A new example is a script in `examples/`, run from the repository root with
  the `examples` environment. Benchmarks compare the result with an analytical solution or
  with another code, and cite the source of the reference solution.
- **Documentation.** Exported functions and types have docstrings, including the units of
  physical quantities (SI units are used throughout).
- **Style.** The code is formatted with JuliaFormatter (Blue style, 120 characters).

## Releases

EKRem follows [Semantic Versioning](https://semver.org/). A release is made by the
maintainer:

1. the version in `Project.toml` is increased (`make bump-patch`, `make bump-minor` or
   `make bump-major`), and the section `[Unreleased]` of the changelog is renamed after the
   new version;
2. a Git tag (`vX.Y.Z`) and a GitHub release with the changelog entry are created.

## Support and governance

EKRem is developed and maintained by Matteo Masi. Issues and pull requests are reviewed on a
best-effort basis, usually within two weeks. Decisions on the scope and the design of the
package are taken by the maintainer, after discussion in the corresponding issue.
Contributors whose work is merged are acknowledged in the release notes.
