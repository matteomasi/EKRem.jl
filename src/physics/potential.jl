# Formulations of the equation for the electric potential.

"Formulation of the equation for the electric potential."
abstract type PotentialModel end

@doc raw"""
    Poisson()

The electric potential ``\phi`` is described by the Poisson equation, written per unit bulk
volume of the porous medium:

```math
-\nabla \cdot (\varepsilon_\mathrm{eff} \nabla \phi) = n F \sum_i z_i c_i,
\qquad \varepsilon_\mathrm{eff} = n \varepsilon_0 \varepsilon_r,
```

where ``c_i`` are the pore-water concentrations. Where the Debye length is much smaller
than the grid spacing, the solution is electroneutral except in boundary layers, and it
does not depend on ``\varepsilon_r``.
"""
struct Poisson <: PotentialModel end

@doc raw"""
    ChargeConservation()

The electric potential ``\phi`` is determined by the conservation of electric charge in an
electroneutral solution:

```math
\nabla \cdot \mathbf{i} = 0, \qquad \mathbf{i} = F \sum_i z_i \mathbf{J}_i,
```

where ``\mathbf{J}_i`` are the Nernst–Planck fluxes (per unit bulk area) and
``\mathbf{i}`` is the current density. The current density includes the electromigration current,
proportional to the conductivity ``\sigma = F^2 / (R T) \sum_i z_i^2 D_i^* c_i``, and the
diffusion current.

This is the limit of the [`Poisson`](@ref) formulation for vanishing Debye length. The
solution does not depend on the permittivity. The initial and boundary compositions must
be electroneutral; electroneutrality is then maintained, because the sum of the transport
equations weighted by ``z_i`` gives
``n \, \partial (\sum_i z_i c_i) / \partial t = -\nabla \cdot \mathbf{i} / F = 0``.
This formulation also allows a fixed current density to be applied at one boundary.
"""
struct ChargeConservation <: PotentialModel end

@doc raw"""
    FixedField()

The electric potential is not coupled to the concentrations. It is the solution of the
Laplace equation ``\nabla^2 \phi = 0`` with the potentials of the boundaries, that is, the
potential in a medium of uniform conductivity. In one dimension, the electric field is
uniform, ``E = (\phi_\mathrm{left} - \phi_\mathrm{right}) / L``, and constant in time.

This formulation is used for comparison with analytical solutions and with models that
assume a constant electric field. The ion transport then does not affect the current, and
the current density is not conserved.
"""
struct FixedField <: PotentialModel end
