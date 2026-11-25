# Geochemistry

## PhreeqcRM library

```@docs
set_phreeqcrm_library
phreeqcrm_library
PhreeqcRMError
```

## Reactions in a set of cells

```@docs
Chemistry
aqueous_species
species_concentrations
set_species_concentrations!
component_concentrations
react!
selected_output(::Chemistry)
```

## PhreeqcRM functions

```@docs
EKRem.LibPhreeqcRM
```

The functions are described in the
[PhreeqcRM C reference](https://phreeqc-dev.github.io/phreeqcrm/RM__interface__C_8h.html).
The following functions are available:

```@eval
using EKRem, Markdown
names_rm = sort(filter(n -> startswith(string(n), "RM_"), names(EKRem.LibPhreeqcRM; all=true)))
Markdown.parse(join(("`$n`" for n in names_rm), ", ") * ".")
```
