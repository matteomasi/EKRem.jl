# Formulations of the equation for the electric potential.

"Formulation of the equation for the electric potential."
abstract type PotentialModel end

"""
    Poisson()

The electric potential φ is described by the Poisson equation, written per unit bulk
volume of the porous medium:

    −∇·(ε_eff ∇φ) = n F Σᵢ zᵢ cᵢ,    ε_eff = n ε₀ εr,

where cᵢ are the pore-water concentrations. Where the Debye length is much smaller than
the grid spacing, the solution is electroneutral except in boundary layers, and it does
not depend on εr.
"""
struct Poisson <: PotentialModel end

"""
    ChargeConservation()

The electric potential φ is determined by the conservation of electric charge in an
electroneutral solution:

    ∇·i = 0,    i = F Σᵢ zᵢ Jᵢ,

where Jᵢ are the Nernst–Planck fluxes (per unit bulk area) and i is the current density.
The current density includes the migration current, proportional to the conductivity
σ = F²/(R T) Σᵢ zᵢ² Dᵢ* cᵢ, and the diffusion current.

This is the limit of the [`Poisson`](@ref) formulation for vanishing Debye length. The
solution does not depend on the permittivity. The initial and boundary compositions must
be electroneutral; electroneutrality is then maintained, because the sum of the transport
equations weighted by zᵢ gives n ∂(Σᵢ zᵢ cᵢ)/∂t = −∇·i / F = 0. This formulation also
allows a fixed current density to be applied at one boundary.
"""
struct ChargeConservation <: PotentialModel end

"""
    FixedField()

The electric potential is not coupled to the concentrations. It is the solution of the
Laplace equation ∇²φ = 0 with the potentials of the boundaries, that is, the potential in
a medium of uniform conductivity. In one dimension, the electric field is uniform,
E = (φ_left − φ_right) / L, and constant in time.

This formulation is used for comparison with analytical solutions and with models that
assume a constant electric field. The ion transport then does not affect the current, and
the current density is not conserved.
"""
struct FixedField <: PotentialModel end
