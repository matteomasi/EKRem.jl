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
