# Private functions

Functions and types that are not exported. They are described here for development and may
change without notice.

## Transport problem

```@docs
SimEK.TransportProblem
SimEK.NPPData
SimEK.solve_transport
SimEK.solve_potential
SimEK.solve_potential!
SimEK.ode_function
SimEK.residual!
SimEK.update_darcy!
SimEK.electrode_production!
SimEK.SOLVERS
SimEK.PROGRESS_INTERVAL
```

## Species and medium

```@docs
SimEK.concentration_vector
SimEK.charge_imbalance
SimEK.effective_diffusion
SimEK.effective_permittivity
```

## Boundaries and coupling

```@docs
SimEK.fixed_potential
SimEK.Coupling
SimEK.transported_species
SimEK.react_nodes!
```
