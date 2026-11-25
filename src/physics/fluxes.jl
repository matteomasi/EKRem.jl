# Numerical fluxes. All functions act on a single edge between two control volumes and are
# therefore independent of the spatial dimension.

@doc raw"""
    sg_flux(D, z, c_k, c_l, Δψ, qh=0)

Scharfetter–Gummel approximation of the Nernst–Planck flux from node ``k`` to node ``l``,
multiplied by the edge length ``h_{kl}``:

```math
h_{kl} J_{kl} = D \left[ B(x)\, c_k - B(-x)\, c_l \right],
\qquad x = z \Delta\psi - \frac{q h_{kl}}{D},
\qquad B(x) = \frac{x}{e^x - 1},
```

where `D` is the (effective) diffusion coefficient, `z` the charge number,
``\Delta\psi = F (\phi_l - \phi_k) / (R T)`` the dimensionless potential difference along
the edge and `qh` the product of the water flux ``q`` (Darcy flux from ``k`` to ``l``,
m s⁻¹) and ``h_{kl}``. The expression reduces to central differences for small ``|x|`` and to
upwinding for large ``|x|``; for ``D = 0`` the flux is the upwind advective flux. This
scaling (flux × ``h``) is the convention of VoronoiFVM.jl flux callbacks.

The two Bernoulli coefficients are evaluated separately. `VoronoiFVM.fbernoulli_pm`
obtains ``B(-x)`` as ``x + B(x)``, which loses relative accuracy in the smaller coefficient
for large ``|x|``; the resulting flux error is of the order of the machine precision times
the concentration ratio across the edge, which is not negligible for trace species.
"""
@inline function sg_flux(D, z, c_k, c_l, Δψ, qh=0.0)
    iszero(D) && return max(qh, zero(qh)) * c_k + min(qh, zero(qh)) * c_l
    x = z * Δψ - qh / D
    return D * (fbernoulli(x) * c_k - fbernoulli(-x) * c_l)
end
