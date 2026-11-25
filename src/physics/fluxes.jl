# Numerical fluxes. All functions act on a single edge between two control volumes and are
# therefore independent of the spatial dimension.

@doc raw"""
    sg_flux(D, z, c_k, c_l, Δψ)

Scharfetter–Gummel approximation of the Nernst–Planck flux from node ``k`` to node ``l``,
multiplied by the edge length ``h_{kl}``:

```math
h_{kl} J_{kl} = D \left[ B(z \Delta\psi)\, c_k - B(-z \Delta\psi)\, c_l \right],
\qquad B(x) = \frac{x}{e^x - 1},
```

where `D` is the (effective) diffusion coefficient, `z` the charge number and
``\Delta\psi = F (\phi_l - \phi_k) / (R T)`` the dimensionless potential difference along
the edge. The expression reduces to central differences for small ``|z \Delta\psi|`` and to
upwinding for large ``|z \Delta\psi|``. This scaling (flux × ``h``) is the convention of
VoronoiFVM.jl flux callbacks.

The two Bernoulli coefficients are evaluated separately. `VoronoiFVM.fbernoulli_pm`
obtains ``B(-x)`` as ``x + B(x)``, which loses relative accuracy in the smaller coefficient
for large ``|x|``; the resulting flux error is of the order of the machine precision times
the concentration ratio across the edge, which is not negligible for trace species.
"""
@inline function sg_flux(D, z, c_k, c_l, Δψ)
    x = z * Δψ
    return D * (fbernoulli(x) * c_k - fbernoulli(-x) * c_l)
end
