# Numerical fluxes. All functions act on a single edge between two control volumes and are
# therefore independent of the spatial dimension.

"""
    sg_flux(D, z, c_k, c_l, Δψ)

Scharfetter–Gummel approximation of the Nernst–Planck flux from node k to node l,
multiplied by the edge length h_kl:

    h J_kl = D [B(z Δψ) c_k − B(−z Δψ) c_l],    B(x) = x / (eˣ − 1),

where `D` is the (effective) diffusion coefficient, `z` the charge number and
Δψ = F (φ_l − φ_k) / (R T) the dimensionless potential difference along the edge.
The expression reduces to central differences for small |z Δψ| and to upwinding for
large |z Δψ|. This scaling (flux × h) is the convention of VoronoiFVM.jl flux callbacks.

The two Bernoulli coefficients are evaluated separately. `VoronoiFVM.fbernoulli_pm`
obtains B(−x) as x + B(x), which loses relative accuracy in the smaller coefficient for
large |x|; the resulting flux error is of the order of the machine precision times the
concentration ratio across the edge, which is not negligible for trace species.
"""
@inline function sg_flux(D, z, c_k, c_l, Δψ)
    x = z * Δψ
    return D * (fbernoulli(x) * c_k - fbernoulli(-x) * c_l)
end
