# [Geochemical reactions](@id reactions-guide)

Geochemical reactions are included when a PHREEQC `database` and an `input` are passed to
[`Model1D`](@ref). The reactions are computed by PhreeqcRM in every node and coupled with
transport as described in [Coupling with geochemical reactions](@ref coupling). The
PhreeqcRM library must be set first (see [Installation](installation.md)).

## PHREEQC input

The input contains the definitions of the solutions and reactants (`SOLUTION`,
`EQUILIBRIUM_PHASES`, `KINETICS`, `EXCHANGE`, `SURFACE`, `GAS_PHASE`,
`SOLID_SOLUTIONS`), together with any additional `PHASES`, `SOLUTION_SPECIES`, `RATES` and
`SELECTED_OUTPUT` blocks. It is given as text or as the path of a file. The keywords are
described in the [PHREEQC manual](https://www.usgs.gov/software/phreeqc-version-3).

The definitions are assigned by number:

- `initial` is a named tuple with the numbers of the definitions in the column, for example
  `(solution=1, equilibrium_phases=1)`; each value is one number for all nodes or a vector
  with one number per node;
- a reservoir is given by the number of a solution, `Reservoir(; solution=2, potential=0.0)`;
- a compartment is given by a named tuple like `initial`,
  `Compartment(; initial=(solution=3, equilibrium_phases=3), volume, current_density)`.

Each PhreeqcRM cell contains 1 L of pore water. The amounts of minerals and other reactants
are therefore given **per litre of pore water**, not per kilogram of water.

## Transported species

By default, all aqueous species of the PHREEQC system except H₂O are transported, with the
charge numbers and the diffusion coefficients at 25 °C of the database (`-dw`; species
without a diffusion coefficient do not move). `species` selects a subset, as a vector of
names (coefficients from the database) or of [`Species`](@ref) (coefficients given by the
user). The species that are not transported remain in the cells.

Every element of the system must be present in the input, in a solution or a reactant,
for PhreeqcRM to include its species.

## Example

Calcite in a column flushed by acid produced at the anode; the catholyte is kept at pH 3 by
the addition of nitric acid through a fictitious phase `Fix_H+`:

```julia
input = """
PHASES
Fix_H+
    H+ = H+
    log_k 0
SOLUTION 1          # pore water
    pH 7 charge
    Na 10
    Cl 10
EQUILIBRIUM_PHASES 1
    Calcite 0 0.1
SOLUTION 2          # electrolyte of the compartments
    pH 7 charge
    Na 1
    N(5) 1
EQUILIBRIUM_PHASES 2
    Fix_H+ -3 HNO3 10
SELECTED_OUTPUT 1
    -reset false
    -pH true
END
"""

model = Model1D(;
    length=0.02,
    grid=collect(range(0, 0.02; length=21)),
    initial=(solution=1, equilibrium_phases=1),
    left=Compartment(; initial=(solution=2,), volume=1.0e-4, current_density=10.0),
    right=Compartment(; initial=(solution=2, equilibrium_phases=2), volume=1.0e-4, potential=0.0),
    database="path/to/phreeqc.dat",
    input,
    potential=ChargeConservation(),
    porosity=0.4,
    tortuosity=0.5,
    area=1.0e-3,
)
run!(model, [300.0]; coupling_step=60.0)
pH = selected_output(model, 300.0)["pH"]
```

## Coupling step

`coupling_step` (s) is required by [`run!`](@ref) in models with reactions. Each interval
between output times is divided into equal steps not longer than `coupling_step`. The
splitting error decreases in proportion to the step; a suitable value is found by
repeating the simulation with a shorter step and comparing the results, particularly at
reaction fronts.

`threads` (argument of [`Model1D`](@ref)) sets the number of threads of PhreeqcRM; the
reactions in different nodes are computed in parallel.

## Results

In addition to the concentrations of the transported species, the columns of the PHREEQC
`SELECTED_OUTPUT` block (for example pH, saturation indices or the amounts of minerals) are
stored at every output time and returned by [`selected_output`](@ref). Total amounts in the
column and the compartments are returned by [`amount`](@ref).

The PhreeqcRM instance can also be used directly through [`Chemistry`](@ref), for example to
test an input before it is used in a transport model.
