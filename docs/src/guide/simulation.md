# Simulation and results

## Running a model

[`run!`](@ref) integrates the model from its current time to the largest of the given output
times (s) and stores the solution at each of them. Repeated calls continue from the last
state, so a simulation can be extended or split into parts with different settings:

```julia
run!(model, [3600.0, 86400.0])               # output after 1 h and 1 d
run!(model, [7 * 86400.0]; solver=:FBDF)     # continue up to 7 d with FBDF
```

The time step is adapted automatically. `reltol` is the relative error tolerance of the
concentrations, and `abstol` the absolute tolerance relative to the largest concentration
of each species.

## Time integration methods

| `solver` | Method | Suited to |
|:---------|:-------|:----------|
| `:Rodas5P` (default) | Rosenbrock method of order 5 | sharp fronts that move by migration |
| `:FBDF` | backward differentiation formulas of variable order (1 to 5) | long simulations with slowly varying solutions |

`Rodas5P` conserves the amounts of the species to round-off. `FBDF` conserves them only to
the tolerance of its Newton iteration, with relative errors of the order of `reltol`.

## Results

| Function | Result |
|:---------|:-------|
| [`coordinates`](@ref) | node coordinates (m) |
| [`output_times`](@ref) | times of the stored results (s), starting with 0 |
| [`concentration`](@ref) | pore-water concentration of a species at the nodes (mol m⁻³) |
| [`electric_potential`](@ref) | potential at the nodes (V) |
| [`voltage`](@ref) | potential difference between the boundaries (V) |
| [`current_density`](@ref) | current density between neighbouring nodes (A m⁻²) |
| [`amount`](@ref) | amount of a species in the column and the compartments (mol m⁻²) |
| [`selected_output`](@ref) | PHREEQC selected output at the nodes |

Except for `coordinates` and `output_times`, the functions take the model and an output
time:

```julia
x = coordinates(model)
c = concentration(model, "Na+", 86400.0)
φ = electric_potential(model, 86400.0)
```
