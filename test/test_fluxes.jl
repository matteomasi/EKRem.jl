# Unit tests of the Scharfetter–Gummel flux, sg_flux(D, z, c_k, c_l, Δψ) = h J_kl.
#
# Notation: x = z Δψ is the dimensionless potential-energy difference along the edge, and
# B(x) = x / (eˣ − 1) is the Bernoulli function, with B(−x) = B(x) + x.

using Test
using EKRem

const D = 1.0e-9   # m² s⁻¹

@testset "Scharfetter–Gummel flux" begin
    @testset "Diffusion limit" begin
        # Without a potential difference, or for a neutral species, Fick's law is recovered
        @test sg_flux(D, 1, 3.0, 1.0, 0.0) == D * (3.0 - 1.0)
        @test sg_flux(D, 0, 3.0, 1.0, 25.0) == D * (3.0 - 1.0)
    end

    @testset "Central-difference limit (|zΔψ| ≪ 1)" begin
        # B(±x) = 1 ∓ x/2 + O(x²): h J = −D [(c_l − c_k) + z Δψ (c_k + c_l)/2]
        for z in (-2, -1, 1, 2), Δψ in (-1.0e-4, 1.0e-4)
            c_k, c_l = 3.0, 1.0
            central = -D * ((c_l - c_k) + z * Δψ * (c_k + c_l) / 2)
            @test sg_flux(D, z, c_k, c_l, Δψ) ≈ central rtol = 1.0e-8
        end
    end

    @testset "Upwind limit (|zΔψ| ≫ 1)" begin
        # The flux is carried by the upstream concentration only: for x = zΔψ ≫ 1 the ions
        # move from l to k (h J ≈ −D x c_l); for x ≪ −1 they move from k to l (h J ≈ −D x c_k)
        c_k, c_l = 3.0, 1.0
        @test sg_flux(D, 1, c_k, c_l, 50.0) ≈ -D * 50.0 * c_l rtol = 1.0e-12
        @test sg_flux(D, 1, c_k, c_l, -50.0) ≈ D * 50.0 * c_k rtol = 1.0e-12
        @test sg_flux(D, -2, c_k, c_l, 25.0) ≈ D * 50.0 * c_k rtol = 1.0e-12
    end

    @testset "Exactness for a constant field" begin
        # For a uniform field the steady profile between two nodes is
        #     c(s) = A + C exp(−x s),   s = ξ / h ∈ [0, 1],
        # with exact flux h J = −D x A (own derivation). The Scharfetter–Gummel flux
        # reproduces it for any x and any C.
        for z in (-2, -1, 1, 2), Δψ in (-20.0, -3.0, -0.5, 0.5, 3.0, 20.0), (A, C) in ((1.0, 2.0), (2.0, -0.5))
            x = z * Δψ
            c_k = A + C
            c_l = A + C * exp(-x)
            @test sg_flux(D, z, c_k, c_l, Δψ) ≈ -D * x * A rtol = 1.0e-10
        end
    end

    @testset "Advection" begin
        # With a water flux q, the steady profile between two nodes for a uniform field is
        # c(s) = A + C exp(−x s), with x = zΔψ − q h / D, and h J = −D x A
        for z in (0, 1, -2),
            Δψ in (-3.0, 0.0, 0.5), qh in (-20.0, -1.0, 2.0, 30.0) .* D,
            (A, C) in ((1.0, 2.0), (2.0, -0.5))

            x = z * Δψ - qh / D
            @test sg_flux(D, z, A + C, A + C * exp(-x), Δψ, qh) ≈ -D * x * A rtol = 1.0e-10
        end
        # Without diffusion, the flux is the upwind advective flux
        @test sg_flux(0.0, 1, 3.0, 1.0, 5.0, 2.0) == 2.0 * 3.0
        @test sg_flux(0.0, 1, 3.0, 1.0, 5.0, -2.0) == -2.0 * 1.0
        @test sg_flux(0.0, 1, 3.0, 1.0, 5.0) == 0
    end

    @testset "Zero flux at Boltzmann equilibrium" begin
        # c_l / c_k = exp(−zΔψ) is the equilibrium distribution (special case A = 0 above).
        # Only |zΔψ| ≤ 50 is tested: beyond, VoronoiFVM sets the smaller Bernoulli coefficient
        # to zero, which corresponds to concentration ratios above 5e21 across one edge.
        for z in (-1, 1, 3), Δψ in (-15.0, -1.0, 0.1, 1.0, 15.0)
            c_k = 2.0
            c_l = c_k * exp(-z * Δψ)
            @test abs(sg_flux(D, z, c_k, c_l, Δψ)) <= 1.0e-12 * D * c_k * max(1.0, abs(z * Δψ))
        end
    end

    @testset "Conservation and symmetry" begin
        for z in (-2, 1), Δψ in (-7.0, 0.3, 7.0)
            c_k, c_l = 3.0, 1.0
            # The flux from l to k is the negative of the flux from k to l
            @test sg_flux(D, z, c_l, c_k, -Δψ) ≈ -sg_flux(D, z, c_k, c_l, Δψ)
            # Reversing the charge is equivalent to reversing the potential difference
            @test sg_flux(D, -z, c_k, c_l, Δψ) ≈ sg_flux(D, z, c_k, c_l, -Δψ)
            # The flux is linear in D
            @test sg_flux(2D, z, c_k, c_l, Δψ) ≈ 2 * sg_flux(D, z, c_k, c_l, Δψ)
        end
    end

    @testset "Monotonicity" begin
        # The flux is non-decreasing in c_k and non-increasing in c_l for every field
        # strength. This property (non-negative Bernoulli coefficients) gives the
        # non-negativity of the concentrations with implicit time integration. For |x| > 50
        # the downstream coefficient is set to zero by VoronoiFVM (B(50) ≈ 1e-20).
        for x in (-700.0, -50.0, -1.0, -1.0e-8, 0.0, 1.0e-8, 1.0, 50.0, 700.0)
            @test sg_flux(D, 1, 1.0, 0.0, x) >= 0
            @test sg_flux(D, 1, 0.0, 1.0, x) <= 0
        end
        for x in (-40.0, -1.0, 0.0, 1.0, 40.0)
            @test sg_flux(D, 1, 1.0, 0.0, x) > 0
            @test sg_flux(D, 1, 0.0, 1.0, x) < 0
        end
    end

    @testset "Accuracy for large concentration ratios" begin
        # Concentration ratio e^|x| across one edge (5e8 for |x| = 20, 2e17 for |x| = 40),
        # slightly away from Boltzmann equilibrium, so that the two terms of the flux nearly
        # cancel. The reference is the same expression evaluated in BigFloat arithmetic.
        bernoulli_big(x) = x / expm1(x)
        for x in (-40.0, -20.0, 20.0, 40.0)
            c_k = 1.0
            c_l = c_k * exp(-x) * (1 + 1.0e-6)
            reference = D * (bernoulli_big(big(x)) * big(c_k) - bernoulli_big(big(-x)) * big(c_l))
            @test sg_flux(D, 1, c_k, c_l, x) ≈ Float64(reference) rtol = 1.0e-8
        end
    end
end
