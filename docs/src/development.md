# Development set-up

The common development tasks are run with `make` (GNU Make) from the repository root;
`make help` lists them:

```console
make install-dev    # instantiate the environments and install the development tools
make test           # run the unit tests (directly in the package environment)
make test-pkg       # run the unit tests with Pkg.test, in an isolated test environment
make examples       # run all examples (figures in examples/output)
make docs           # build the documentation (site in docs/build/1)
make docs-serve     # build the documentation and serve it locally with live reload
make lint           # check package quality with Aqua.jl
make format         # format the code with JuliaFormatter (settings in .JuliaFormatter.toml)
```

Without `make`, the environment is instantiated and the tests are run from the Julia
package manager:

```julia
pkg> activate .
pkg> instantiate
pkg> test
```

## Tests

The tests in `test/` are fast unit tests that check results following from physics and
chemistry: conservation of mass and charge, equilibrium, stoichiometry and analytical
solutions. They do not compare with values that may change between versions of a database
or a library. Complete simulations belong in `examples/`.

Without a PhreeqcRM library, the tests of the geochemical reactions are skipped. By
default, the tests read the database `../database/phreeqc.dat` relative to the library;
another file is set with the environment variable `EKREM_PHREEQC_DATABASE`.

## Documentation

The documentation is built with [Documenter.jl](https://documenter.juliadocs.org/) and
[DocumenterVitepress.jl](https://github.com/LuxDL/DocumenterVitepress.jl). The pages are in
`docs/src`; the example pages are generated from the scripts in `examples/` by
[Literate.jl](https://fredrikekre.github.io/Literate.jl/), and the examples are run during the
build. In the scripts, lines ending with `#src` are omitted from the pages, and lines
starting with `#md` appear only in the pages.

Equations are written in LaTeX: inline as ``` ``\sigma = F^2 / (R T)`` ```, and in display
mode in a `math` code block. Docstrings that contain backslashes are written as
`@doc raw"""..."""`.

The first build installs the VitePress packages with npm, which requires internet access.
