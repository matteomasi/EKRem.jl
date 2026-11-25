# Private functions

Functions and types that are not exported. They are described here for development and may
change without notice.

## Transport problem

```@docs
EKRem.TransportProblem
EKRem.NPPData
EKRem.solve_transport
EKRem.solve_potential
EKRem.solve_potential!
EKRem.ode_function
EKRem.residual!
EKRem.SOLVERS
EKRem.PROGRESS_INTERVAL
```

## Species and medium

```@docs
EKRem.concentration_vector
EKRem.charge_imbalance
EKRem.effective_diffusion
EKRem.effective_permittivity
```

## Boundaries and coupling

```@docs
EKRem.fixed_potential
EKRem.Coupling
EKRem.transported_species
EKRem.react_nodes!
```
