# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/fite.py

using StaticArrays

# ============================================================================
# Abstract Interpolator Interface
# ============================================================================

"""
    AbstractInterpolator

Abstract type for band interpolation backends. It is the interpolator axis of
[`solve_transport`](@ref) and of `run_integrate(interp, sys; sampling)`.

Concrete implementations:

  - [`FourierInterpolator`](@ref): BoltzTraP-style Fourier interpolation
    (eigenvalues only)
  - [`WannierInterpolator`](@ref): Wannier.jl-based interpolation (provided by
    the Wannier.jl weakdep extension)

Common API:

  - `interpolate_bands(interp, kpoints)` -> eigenvalues
  - [`interpolate_velocities`](@ref)`(interp, kpoints)` -> group velocities
"""
abstract type AbstractInterpolator end

#=
    interpolate_bands(interp::AbstractInterpolator, kpoints)

Interpolate band energies at given k-points.

# Arguments
- `interp`: Interpolator instance
- `kpoints`: Matrix of k-points (nk × 3) in fractional coordinates

# Returns
- `ebands`: Band energies (nbands × nk)
=#
function interpolate_bands end

"""
    interpolate_velocities(interp::AbstractInterpolator, kpoints)

Interpolate band group velocities at given k-points.

# Arguments

  - `interp`: an [`AbstractInterpolator`](@ref)
  - `kpoints`: Matrix of k-points (nk × 3) in fractional coordinates

# Returns

  - `vbands`: Group velocities (3 × nbands × nk), in Hartree·Bohr

These are the band-diagonal group velocities. The degeneracy-aware,
gauge-invariant velocity tensor used inside the transport pipeline is built
separately from the full velocity matrices.
"""
function interpolate_velocities end

#=
    velocity_matrices(interp, kpoints)

Interpolate the full band-gauge covariant velocity matrices (n_wann × n_wann × 3 ×
nk) at given k-points. Unlike [`interpolate_velocities`] (diagonal group velocity
only), this exposes the intra-degenerate-manifold off-diagonal velocity matrix
elements needed for a gauge-invariant transport tensor at degenerate k-points.

Internal extension point: only the Wannier path implements it (in the Wannier
extension); it is not part of the public `AbstractInterpolator` contract.
=#
function velocity_matrices end

#=
    interpolate(interp::AbstractInterpolator, kpoints)

Interpolate both band energies and velocities at given k-points.

# Returns
- `ebands`: Band energies (nbands × nk)
- `vbands`: Group velocities (3 × nbands × nk)
=#
function interpolate(interp::AbstractInterpolator, kpoints)
    ebands = interpolate_bands(interp, kpoints)
    vbands = interpolate_velocities(interp, kpoints)
    return ebands, vbands
end

# Error dispatch for non-AbstractInterpolator types
function interpolate_bands(interp, kpoints)
    throw(
        ArgumentError(
            "interpolate_bands expects a concrete subtype of AbstractInterpolator, got $(typeof(interp)).",
        ),
    )
end

function interpolate_velocities(interp, kpoints)
    throw(
        ArgumentError(
            "interpolate_velocities expects a concrete subtype of AbstractInterpolator, got $(typeof(interp)).",
        ),
    )
end

function interpolate(interp, kpoints)
    throw(
        ArgumentError(
            "interpolate expects a concrete subtype of AbstractInterpolator, got $(typeof(interp)).",
        ),
    )
end

# Abstract-type fallback for AbstractInterpolator subtypes that lack a
# concrete implementation. Without these, dispatch falls to the ::Any
# fallback above and the error message states "expects a concrete subtype"
# even when the caller passed a proper subtype, which is misleading. (Issue #54.)
#
# Note: `interpolate(::AbstractInterpolator, kpoints)` already exists above
# (line 64) as the default delegate that calls `interpolate_bands` and
# `interpolate_velocities`. For a subtype-without-impl, that delegate
# cascades through the new abstract-type methods below, so a separate
# abstract-type method for `interpolate` is not needed.
function interpolate_bands(interp::AbstractInterpolator, kpoints)
    error(
        "interpolate_bands: not implemented for $(typeof(interp)). " *
        "Subtypes of AbstractInterpolator must define their own method.",
    )
end

function interpolate_velocities(interp::AbstractInterpolator, kpoints)
    error(
        "interpolate_velocities: not implemented for $(typeof(interp)). " *
        "Subtypes of AbstractInterpolator must define their own method.",
    )
end

# ============================================================================
# Convenience Functions (Default to FourierInterpolator)
# ============================================================================

#=
    interpolate(kpoints, energies, equivalences, lattvec, new_kpoints)

Convenience function that creates a FourierInterpolator and interpolates at new k-points.

# Arguments
- `kpoints`: Original k-points (nk × 3)
- `energies`: Band energies at original k-points (nbands × nk)
- `equivalences`: Vector of equivalence class matrices
- `lattvec`: 3×3 lattice vectors
- `new_kpoints`: K-points to interpolate at (nk_new × 3)

# Returns
- `ebands`: Interpolated band energies (nbands × nk_new)
- `vbands`: Group velocities (3 × nbands × nk_new)
=#
function interpolate(
    kpoints::AbstractMatrix,
    energies::AbstractMatrix,
    equivalences::Vector,
    lattvec::AbstractMatrix,
    new_kpoints::AbstractMatrix,
)
    interp = FourierInterpolator(kpoints, energies, equivalences, lattvec)
    return interpolate(interp, new_kpoints)
end

#=
    interpolate_bands(kpoints, energies, equivalences, lattvec, new_kpoints)

Convenience function that creates a FourierInterpolator and interpolates band energies.

# Returns
- `ebands`: Interpolated band energies (nbands × nk_new)
=#
function interpolate_bands(
    kpoints::AbstractMatrix,
    energies::AbstractMatrix,
    equivalences::Vector,
    lattvec::AbstractMatrix,
    new_kpoints::AbstractMatrix,
)
    interp = FourierInterpolator(kpoints, energies, equivalences, lattvec)
    return interpolate_bands(interp, new_kpoints)
end

#=
    interpolate_velocities(kpoints, energies, equivalences, lattvec, new_kpoints)

Convenience function that creates a FourierInterpolator and interpolates group velocities.

# Returns
- `vbands`: Group velocities (3 × nbands × nk_new)
=#
function interpolate_velocities(
    kpoints::AbstractMatrix,
    energies::AbstractMatrix,
    equivalences::Vector,
    lattvec::AbstractMatrix,
    new_kpoints::AbstractMatrix,
)
    interp = FourierInterpolator(kpoints, energies, equivalences, lattvec)
    return interpolate_velocities(interp, new_kpoints)
end

# ============================================================================
# Fourier Interpolator (BoltzTraP method)
# ============================================================================

"""
    FourierInterpolator <: AbstractInterpolator

BoltzTraP-style Fourier interpolation using star functions.

# Fields

  - `coeffs`: Fourier coefficients (nbands × neq)
  - `equivalences`: Vector of equivalence class matrices
  - `lattvec`: 3×3 lattice vectors (columns), in Bohr

# Constructors

    FourierInterpolator(kpoints, energies, equivalences, lattvec)
    FourierInterpolator(interp::InterpolationResult)

The first fits band energies to a Fourier series. The second reuses the
coefficients of an existing [`InterpolationResult`](@ref).

The star functions fix the smallest grid that represents them exactly, so this
interpolator derives its own integration grid and ignores an explicit
[`UniformMesh`](@ref) size.
"""
struct FourierInterpolator <: AbstractInterpolator
    coeffs::Matrix{ComplexF64}
    equivalences::Vector{<:AbstractMatrix{<:Integer}}
    lattvec::SMatrix{3,3,Float64}
end

#=
    FourierInterpolator(kpoints, energies, equivalences, lattvec)

Construct a FourierInterpolator by fitting band energies.

# Arguments
- `kpoints`: Matrix of k-points (nk × 3)
- `energies`: Matrix of band energies (nbands × nk)
- `equivalences`: Vector of equivalence class matrices
- `lattvec`: 3×3 lattice vectors
=#
function FourierInterpolator(kpoints, energies, equivalences, lattvec)
    coeffs = fitde3D(kpoints, energies, equivalences, lattvec)
    lattvec_static = SMatrix{3,3,Float64}(lattvec)
    return FourierInterpolator(coeffs, equivalences, lattvec_static)
end

# Number of bands
nbands(interp::FourierInterpolator) = size(interp.coeffs, 1)

# Implement abstract interface
function interpolate_bands(interp::FourierInterpolator, kpoints)
    ebands, _ =
        getBands(kpoints, interp.equivalences, Matrix(interp.lattvec), interp.coeffs)
    return ebands
end

function interpolate_velocities(interp::FourierInterpolator, kpoints)
    _, vbands =
        getBands(kpoints, interp.equivalences, Matrix(interp.lattvec), interp.coeffs)
    return vbands
end

function interpolate(interp::FourierInterpolator, kpoints)
    return getBands(kpoints, interp.equivalences, Matrix(interp.lattvec), interp.coeffs)
end

# ============================================================================
# Low-level Functions (used internally)
# ============================================================================

#=
    compute_phase_factors(kpoints, equivalences)

Compute phase factors exp(2πi k·R) for each k-point and equivalence class.

# Arguments
- `kpoints`: Matrix of k-points (nk × 3)
- `equivalences`: Vector of equivalence class matrices

# Returns
- `phase`: Complex matrix (nk × neq)
=#
function compute_phase_factors(kpoints, equivalences)
    nk = size(kpoints, 1)
    neq = length(equivalences)
    phase = zeros(ComplexF64, nk, neq)

    tpii = 2π * im

    for (j, equiv) in enumerate(equivalences)
        nstar = size(equiv, 1)
        # Vectorized: matrix multiply (nk, 3) * (3, nstar) -> (nk, nstar)
        dot_products = kpoints * equiv'
        # exp and sum over stars (dims=2)
        phase_sum = sum(exp.(tpii .* dot_products), dims = 2)
        phase[:, j] = vec(phase_sum) ./ nstar
    end

    return phase
end

#=
    compute_regularization_weights(equivalences, lattvec; C1=0.75, C2=0.75)

Compute regularization weights for Fourier fitting.

Uses the form: ρ(x) = 1 / ((1 - C1*x²)² + C2*x⁶)
where x = |R| / |R₁| (normalized by first non-zero R).
=#
function compute_regularization_weights(equivalences, lattvec; C1 = 0.75, C2 = 0.75)
    # Get representative vectors and their norms
    Rvecs = [equiv[1, :] for equiv in equivalences]
    norms = [norm(lattvec' * R) for R in Rvecs]

    # Normalize by first non-zero norm
    norm1 = norms[2]  # First non-zero (index 2, since index 1 is origin)
    X2 = (norms ./ norm1) .^ 2

    rhoi = @. 1.0 / ((1.0 - C1 * X2)^2 + C2 * X2^3)
    rhoi[1] = 0.0  # Origin has zero weight

    return rhoi
end

#=
    fitde3D(kpoints, energies, equivalences, lattvec)

Fit band energies to Fourier series using regularized least squares.

Port of BoltzTraP2/fite.py fitde3D function.

# Arguments
- `kpoints`: Matrix of k-points (nk × 3)
- `energies`: Matrix of band energies (nbands × nk)
- `equivalences`: Vector of equivalence class matrices
- `lattvec`: 3×3 lattice vectors

# Returns
- `coeffs`: Fourier coefficients (nbands × neq)
=#
function fitde3D(kpoints, energies, equivalences, lattvec)
    nk = size(kpoints, 1)
    nbands = size(energies, 1)
    neq = length(equivalences)

    # Compute phase factors and regularization weights
    phase = compute_phase_factors(kpoints, equivalences)
    rhoi = compute_regularization_weights(equivalences, lattvec)

    # Energy differences (relative to last k-point)
    # De has shape (nbands, nk-1) in Julia, corresponds to De.T in Python
    De = energies[:, 1:(end-1)] .- energies[:, end:end]

    # Phase differences: phaseR (nk-1, neq)
    phaseR = phase[1:(end-1), :] .- phase[end:end, :]

    # Build regularized matrix (matches Python exactly)
    # Python: Hmat = (phaseR[:, 1:] @ (phaseR[:, 1:] * rhoi[1:]).conj().T).real
    # phaseR[:, 2:end] is (nk-1, neq-1)
    # Element-wise multiply with rhoi (broadcast along rows)
    weighted_phaseR = phaseR[:, 2:end] .* rhoi[2:end]'  # (nk-1, neq-1)

    # Hmat = phaseR[:, 2:end] @ weighted_phaseR.conj().T = (nk-1, nk-1)
    Hmat = real(phaseR[:, 2:end] * conj(weighted_phaseR)')

    # Solve least squares: Hmat @ rlambda = De.T
    # Python: rlambda = sp.linalg.lstsq(Hmat, De)[0]  (SVD, rank-deficient OK)
    # De.T in Python is (nk-1, nbands), our De is (nbands, nk-1)
    # So we need De' for the solve.
    # Julia's `\` on a square matrix uses an LU factorization, which raises
    # SingularException for a rank-deficient Hmat (high-symmetry small cells
    # where the source k-mesh is sparse relative to the equivalence stars).
    # `pinv` (Moore-Penrose, SVD-based) reproduces scipy.linalg.lstsq's
    # minimum-norm least-squares solution and is identical to `\` for the
    # full-rank case up to round-off.
    rlambda = pinv(Hmat) * De'  # (nk-1, nbands)

    # Recover coefficients
    # Python: coeffs = rhoi * (rlambda.T @ phaseR)
    # rlambda.T is (nbands, nk-1), phaseR is (nk-1, neq)
    # rlambda.T @ phaseR = (nbands, neq)
    coeffs = rhoi' .* (rlambda' * phaseR)  # (nbands, neq)

    # First coefficient: E₀ = E_ref - Σ c_R exp(2πi k_ref · R)
    # Python: coeffs[:, 0] = ene.T[-1] - coeffs[:, 1:] @ phase[-1, 1:]
    coeffs[:, 1] = energies[:, end] - coeffs[:, 2:end] * phase[end, 2:end]

    return coeffs
end

#=
    getBands(kpoints, equivalences, lattvec, coeffs)

Sample the energy bands at particular k-points from interpolation coefficients.

Port of BoltzTraP2/fite.py getBands function.

# Arguments
- `kpoints`: Matrix of k-points (nk × 3)
- `equivalences`: Vector of equivalence class matrices
- `lattvec`: 3×3 lattice vectors
- `coeffs`: Interpolation coefficients (nbands × neq)

# Returns
- `ebands`: Reconstructed band energies (nbands × nk)
- `vbands`: Group velocities (3 × nbands × nk)
=#
function getBands(kpoints, equivalences, lattvec, coeffs)
    nk = size(kpoints, 1)
    nbands = size(coeffs, 1)
    neq = length(equivalences)

    tpii = 2π * im

    # Compute phase factors for each equivalence class
    phase = zeros(ComplexF64, neq, nk)  # (neq, nk) - transposed vs compute_phase_factors
    phaseR = zeros(ComplexF64, neq, nk, 3)  # (neq, nk, 3)

    for (j, equiv) in enumerate(equivalences)
        ns = size(equiv, 1)

        # Vectorized: (nk, 3) * (3, ns) -> (nk, ns)
        dot_products = kpoints * equiv'
        phase0 = exp.(tpii .* dot_products)

        # Sum over stars, normalized
        phase[j, :] = vec(sum(phase0, dims = 2)) ./ ns

        # Velocity phase factor: 1j * lattvec @ equiv.T
        # vv = (3, ns)
        vv = im * lattvec * equiv'

        # phaseR[j, :, :] = (phase0 @ vv.T) / nstar
        # (nk, ns) * (ns, 3) -> (nk, 3)
        phaseR[j, :, :] = (phase0 * vv') ./ ns
    end

    # Compute bands: egrid = coeffs @ phase
    # coeffs is (nbands, neq), phase is (neq, nk)
    egrid = real(coeffs * phase)  # (nbands, nk)

    # Compute velocities: vgrid = phaseR.T @ coeffs.T
    # phaseR is (neq, nk, 3), coeffs is (nbands, neq)
    # We want (3, nbands, nk)
    # Note: Negate to match Python sign convention
    vgrid = zeros(3, nbands, nk)
    for idir = 1:3
        # phaseR[:, :, idir] is (neq, nk)
        vgrid[idir, :, :] = -real(coeffs * phaseR[:, :, idir])
    end

    return egrid, vgrid
end
