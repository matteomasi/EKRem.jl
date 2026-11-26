# Theory

This page summarises the equations solved by EKRem and their numerical solution. The
porous medium is water-saturated, isothermal and rigid. The only water flow is
electroosmotic; hydraulic gradients are not considered.

## [Transport](@id transport)

The unknowns are the pore-water concentrations ``c_i`` (mol m⁻³) of the ``N`` transported
aqueous species and the electric potential ``\phi`` (V). The mass balance of each species,
per unit bulk volume of the porous medium, is

```math
n \frac{\partial c_i}{\partial t} + \nabla \cdot \mathbf{J}_i = 0, \qquad i = 1, \dots, N,
```

with the Nernst–Planck flux (per unit bulk area)

```math
\mathbf{J}_i = -D_i^* \left( \nabla c_i + z_i c_i \nabla \psi \right) + \mathbf{q} c_i,
\qquad \psi = \frac{F \phi}{R T},
\qquad D_i^* = n \tau D_i,
```

where ``n`` is the porosity, ``\tau`` the tortuosity factor (``0 < \tau \le 1``), ``D_i``
the diffusion coefficient in free water, ``z_i`` the charge number, ``F`` the Faraday
constant, ``R`` the gas constant and ``T`` the temperature. The terms of the flux describe
diffusion, electromigration and advection by the electroosmotic water flux ``\mathbf{q}``
(see [Electroosmosis](@ref electroosmosis)). Geochemical reactions are not part of
these equations; they are computed separately (see
[Coupling with geochemical reactions](@ref coupling)).

!!! note "Tortuosity factor"
    The tortuosity factor multiplies the diffusion coefficient. Where the effective diffusion
    coefficient is defined as ``n D_i / \tau'`` with a tortuosity ``\tau' \ge 1``, the
    corresponding factor is ``\tau = 1 / \tau'``.

## [Electric potential](@id electric-potential)

The potential is described by one of three formulations, selected by the `potential`
argument of the model.

**Poisson equation** ([`Poisson`](@ref)). The charge density of the pore water determines
the potential:

```math
-\nabla \cdot (\varepsilon_\mathrm{eff} \nabla \phi) = n F \sum_i z_i c_i,
\qquad \varepsilon_\mathrm{eff} = n \varepsilon_0 \varepsilon_r .
```

Charge separation is confined to layers whose thickness is of the order of the Debye length

```math
\lambda_D = \sqrt{\frac{\varepsilon_0 \varepsilon_r R T}{F^2 \sum_i z_i^2 c_i}},
```

which is of the order of nanometres in natural pore waters. Where ``\lambda_D`` is much
smaller than the grid spacing ([`debye_ratio`](@ref)), the solution is electroneutral
outside thin boundary layers.

**Charge conservation** ([`ChargeConservation`](@ref)). In the limit of vanishing Debye
length, the solution is electroneutral, ``\sum_i z_i c_i = 0``, and the potential follows
from the conservation of charge:

```math
\nabla \cdot \mathbf{i} = 0, \qquad
\mathbf{i} = F \sum_i z_i \mathbf{J}_i
= -\sigma \nabla \phi - F \sum_i z_i D_i^* \nabla c_i,
\qquad \sigma = \frac{F^2}{R T} \sum_i z_i^2 D_i^* c_i ,
```

where ``\mathbf{i}`` is the current density and ``\sigma`` the bulk electric conductivity.
The weighted sum of the mass balances, ``n \, \partial (\sum_i z_i c_i) / \partial t =
-\nabla \cdot \mathbf{i} / F = 0``, shows that electroneutral initial and boundary
compositions remain electroneutral. This formulation is required when a current density is
applied.

**Fixed field** ([`FixedField`](@ref)). The potential is the solution of the Laplace
equation ``\nabla^2 \phi = 0`` with the boundary potentials. It does not depend on the
concentrations; in one dimension, the electric field is uniform and constant in time. This
formulation is used for comparison with analytical solutions that assume a constant field.

## [Boundary conditions](@id boundary-conditions)

**Reservoirs** ([`Reservoir`](@ref)). The concentrations at the boundary are fixed,
``c_i = c_i^\mathrm{res}``, and the composition is electroneutral.

**Electric conditions.** Under voltage control, the potential is fixed at both boundaries.
Under current control, the potential is fixed at one boundary, and a current density
``i_\mathrm{app}`` is applied at the other one:

```math
\mathbf{i} \cdot \boldsymbol{\nu}_\mathrm{in} = i_\mathrm{app},
```

where ``\boldsymbol{\nu}_\mathrm{in}`` is the unit normal pointing into the domain. The
applied current density is positive at the anode, where the current enters the domain, and
negative at the cathode. In one dimension, the current density is then uniform along the
column, and the potential of the boundary is part of the solution.

## [Electroosmosis](@id electroosmosis)

The electroosmotic water flux (Darcy flux, m s⁻¹) is

```math
\mathbf{q} = -k_{eo} \nabla \phi,
```

where ``k_{eo}`` (m² V⁻¹ s⁻¹) is the electroosmotic permeability of the medium
([`PorousMedium`](@ref)). By the Helmholtz–Smoluchowski equation with the porosity and
tortuosity terms introduced by Casagrande, ``k_{eo} = -\varepsilon \zeta n \tau / \mu``,
where ``\zeta`` is the zeta potential of the solid and ``\mu`` the viscosity of water (Masi
et al., 2017, Eq. 6). For the negative zeta potential of most soils, ``k_{eo} > 0`` and water
flows from the anode to the cathode.

In one dimension, the conservation of the water volume requires a uniform flux. Where the
field varies along the column, a pressure gradient develops; with equal hydraulic heads at
both ends, the flux is the mean of ``-k_{eo} \partial\phi/\partial x`` over the column
(Masi et al., 2017, Eq. 8):

```math
q = \frac{k_{eo}}{L} \left[ \phi(0) - \phi(L) \right].
```

The flux is computed from the boundary potentials at the start of each time step and held
fixed over the step; under voltage control it is constant. In an electroneutral solution,
the water flux carries no net charge and does not contribute to the current.

## [Electrode compartments](@id electrode-compartments)

An electrode compartment ([`Compartment`](@ref)) is a well-mixed solution of volume ``V`` in
contact with one end of the column, of cross-sectional area ``A``. It is flushed with
deionised water at the volumetric rate ``Q``, and its concentrations change by the fluxes
exchanged with the column, by the outflow and by water electrolysis:

```math
V \frac{\mathrm{d} c_i}{\mathrm{d} t} = \mp A J_i - (Q + A q_\mathrm{in}) c_i + R_i,
```

where ``J_i`` is the flux at the boundary in the ``+x`` direction (upper sign at the left
boundary, lower sign at the right boundary) and ``q_\mathrm{in} \ge 0`` the electroosmotic
water flux into the compartment, which overflows. The compartment from which water flows
into the column is refilled with deionised water. Electrolysis produces H⁺ at the anode and OH⁻
at the cathode,

```math
\mathrm{2\,H_2O \to O_2 + 4\,H^+ + 4\,e^-}, \qquad
\mathrm{2\,H_2O + 2\,e^- \to H_2 + 2\,OH^-},
```

at the rate ``R = I / F``, where ``I = |i_\mathrm{app}| A`` is the current; the gases are
assumed to escape. Because the rate of electrolysis is set by the current, compartments
require current control. The anode is the boundary where the current enters the column:
the boundary with the applied current density if ``i_\mathrm{app} > 0``, and the other
boundary otherwise.

In the discrete equations, the compartment is represented by the boundary node, with the
storage ``V/A`` per unit area. The control volume of this node also contains half of the
adjacent cell of the column, whose pore volume is usually negligible compared with ``V``.

## [Discretisation](@id discretisation)

The equations are discretised in space by the vertex-centred finite-volume method of
[VoronoiFVM.jl](https://github.com/WIAS-PDELib/VoronoiFVM.jl): each node is surrounded by a
control volume, and the fluxes are evaluated on the edges between neighbouring nodes. The
Nernst–Planck flux from node ``k`` to node ``l`` at distance ``h_{kl}`` is approximated by
the Scharfetter–Gummel expression ([`sg_flux`](@ref))

```math
J_{kl} = \frac{D^*}{h_{kl}} \left[ B(x_{kl})\, c_k - B(-x_{kl})\, c_l \right],
\qquad x_{kl} = z \Delta\psi_{kl} - \frac{q_{kl} h_{kl}}{D^*},
\qquad B(x) = \frac{x}{e^x - 1},
```

with ``\Delta\psi_{kl} = \psi_l - \psi_k`` and ``q_{kl}`` the water flux from ``k`` to ``l``.
The expression is exact for a constant flux, a linear potential and a constant water flux
along the edge. It reduces to central differences where diffusion dominates and to
upwinding where electromigration or advection dominates.
The examples show second-order convergence of the steady-state solutions with respect to the
grid spacing.

The semi-discrete system is a set of ordinary differential equations for the
concentrations, coupled with algebraic equations for the potential (which has no storage
term). It is integrated in time with adaptive steps and error control, by the Rosenbrock
method `Rodas5P` or by the variable-order backward differentiation formulas `FBDF` of
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) (see [`run!`](@ref)). The
potential is excluded from the error control, because its accuracy follows from that of the
concentrations.

## [Coupling with geochemical reactions](@id coupling)

Geochemical reactions are computed by
[PhreeqcRM](https://phreeqc-dev.github.io/phreeqcrm/) in every node and coupled with
transport by the sequential non-iterative approach. The simulated period is divided into
coupling steps ``\Delta t_c``, and over each step:

1. the transport equations are integrated from ``t`` to ``t + \Delta t_c`` without
   reactions;
2. the concentrations of the transported species are passed to PhreeqcRM, which integrates
   the kinetic reactions over ``\Delta t_c`` and establishes equilibrium in every node;
3. the concentrations after the reactions replace those of the transported species.

Minerals, exchangers, surfaces and other reactants remain in the PhreeqcRM cells and are
not transported. The nodes of reservoirs keep their composition; the nodes of electrode
compartments react like the pore water. Small negative concentrations produced by the time
integration are set to zero before they are passed to PhreeqcRM. At the output times, the
potential is computed again for the concentrations after the reactions.

The splitting error decreases in proportion to ``\Delta t_c``. It is largest at sharp
reaction fronts, such as an acid front advancing in a carbonate medium, and the coupling
step is chosen by comparing the results with those of a shorter step.

Each PhreeqcRM cell contains 1 L of pore water; the amounts of minerals and other reactants
in the PHREEQC input are therefore per litre of pore water. Concentrations are exchanged as
concentrations of aqueous species, in mol per m³ of pore water (equal to mmol L⁻¹), so that
the amounts of the components are conserved by the transfer.

## References

Masi, M., Ceccarini, A., Iannelli, R., 2017. Multispecies reactive transport modelling of
electrokinetic remediation of harbour sediments. J. Hazard. Mater. 326, 187–196.
https://doi.org/10.1016/j.jhazmat.2016.12.032
