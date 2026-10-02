# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/fite.py

# Note: lattvec is atomic units (Bohr) throughout; velocities are Ha·Bohr.

#=
    determine_fft_dims(equivalences)

Determine FFT grid dimensions from equivalence classes.
=#
function determine_fft_dims(equivalences)
    allvec = vcat([Matrix(equiv) for equiv in equivalences]...)
    maxabs = maximum(abs.(allvec), dims = 1)
    return Tuple(2 .* vec(maxabs) .+ 1)
end

# ============================================================================
# Optimized getBTPbands implementation using function barrier pattern
# Achieves ~4x speedup by:
# 1. Pre-computing Cartesian coordinates (lattvec * R)
# 2. Filling all 4 grids in a single loop
# 3. Using function barrier for type stability
# ============================================================================

#=
    _flatten_coeffs(bandcoeff, equivalences) -> Vector{ComplexF64}

Flatten band coefficients with normalization for FFT grid filling.
Each coefficient is normalized by the number of stars in its equivalence class.
=#
function _flatten_coeffs(bandcoeff::AbstractVector, equivalences::Vector)
    total_pts = sum(size(e, 1) for e in equivalences)
    flat = Vector{ComplexF64}(undef, total_pts)
    k = 1
    @inbounds for (c, equiv) in zip(bandcoeff, equivalences)
        nstar = size(equiv, 1)
        c_norm = c / nstar
        for _ = 1:nstar
            flat[k] = c_norm
            k += 1
        end
    end
    return flat
end

#=
    _fill_grids!(egrid, vgrid1, vgrid2, vgrid3, flat_coeffs, allvec_t, all_indices)

Fill energy and velocity FFT grids in a single loop.
Uses pre-computed Cartesian coordinates (allvec_t) and indices for efficiency.
=#
function _fill_grids!(
    egrid::Array{ComplexF64,3},
    vgrid1::Array{ComplexF64,3},
    vgrid2::Array{ComplexF64,3},
    vgrid3::Array{ComplexF64,3},
    flat_coeffs::Vector{ComplexF64},
    allvec_t::Matrix{Float64},
    all_indices::Vector{CartesianIndex{3}},
)
    fill!(egrid, 0)
    fill!(vgrid1, 0)
    fill!(vgrid2, 0)
    fill!(vgrid3, 0)

    @inbounds for k in eachindex(all_indices)
        idx = all_indices[k]
        c = flat_coeffs[k]
        egrid[idx] = c
        vgrid1[idx] = im * allvec_t[1, k] * c
        vgrid2[idx] = im * allvec_t[2, k] * c
        vgrid3[idx] = im * allvec_t[3, k] * c
    end
end

#=
    _getBTPbands_kernel!(eband, vvband, coeffs, equivalences, allvec_t, all_indices, dims, plan)

Inner kernel for getBTPbands with concrete types (function barrier pattern).
This function is type-stable and optimized for performance.
=#
function _getBTPbands_kernel!(
    eband::Matrix{Float64},
    vvband::Union{Array{Float64,4},Nothing},
    coeffs::AbstractMatrix,
    equivalences::Vector{<:AbstractMatrix{<:Integer}},
    allvec_t::Matrix{Float64},
    all_indices::Vector{CartesianIndex{3}},
    dims::NTuple{3,Int},
    plan,
)
    nbands = size(coeffs, 1)
    npts = prod(dims)
    scale = prod(dims)

    # Work arrays (reused across bands)
    egrid = zeros(ComplexF64, dims)
    vgrid1 = zeros(ComplexF64, dims)
    vgrid2 = zeros(ComplexF64, dims)
    vgrid3 = zeros(ComplexF64, dims)

    for iband = 1:nbands
        flat_coeffs = _flatten_coeffs(coeffs[iband, :], equivalences)
        _fill_grids!(egrid, vgrid1, vgrid2, vgrid3, flat_coeffs, allvec_t, all_indices)

        eband[iband, :] = scale .* real.(vec(plan * egrid))

        if !isnothing(vvband)
            vb1 = scale .* real.(vec(plan * vgrid1))
            vb2 = scale .* real.(vec(plan * vgrid2))
            vb3 = scale .* real.(vec(plan * vgrid3))

            # Exploit symmetry: vvband[i,j] = vvband[j,i] (6 multiplications instead of 9)
            @inbounds for k = 1:npts
                # Diagonal components (3 multiplications)
                vvband[iband, 1, 1, k] = vb1[k] * vb1[k]
                vvband[iband, 2, 2, k] = vb2[k] * vb2[k]
                vvband[iband, 3, 3, k] = vb3[k] * vb3[k]

                # Off-diagonal components (3 multiplications + 6 assignments)
                v12 = vb1[k] * vb2[k]
                v13 = vb1[k] * vb3[k]
                v23 = vb2[k] * vb3[k]
                vvband[iband, 1, 2, k] = v12
                vvband[iband, 2, 1, k] = v12
                vvband[iband, 1, 3, k] = v13
                vvband[iband, 3, 1, k] = v13
                vvband[iband, 2, 3, k] = v23
                vvband[iband, 3, 2, k] = v23
            end
        end
    end
end

#=
    getBTPbands(coeffs, equivalences, lattvec; compute_velocity=true)

Reconstruct all bands using FFT.

Uses function barrier pattern for performance optimization:
1. Pre-computes Cartesian coordinates (lattvec * R) outside hot loop
2. Fills all 4 grids (energy + 3 velocity) in a single loop
3. Calls type-stable inner kernel

# Returns
- `eband`: Band energies (nbands × npts)
- `vvband`: Velocity outer products (nbands × 3 × 3 × npts) if compute_velocity=true

Note: lattvec is expected in Bohr (atomic units). Velocities are returned in
Ha·Bohr, consistent with the atomic-units convention used throughout BT.jl.
=#
function getBTPbands(coeffs, equivalences, lattvec; compute_velocity = true)
    dims = determine_fft_dims(equivalences)
    nbands = size(coeffs, 1)
    npts = prod(dims)
    d1, d2, d3 = dims

    # Pre-computation (amortized cost)
    plan = plan_ifft(zeros(ComplexF64, dims))

    # Pre-compute Cartesian coordinates: allvec_t[i, k] = (equiv * lattvec)[k, i]
    all_equiv_mat = vcat([Matrix(equiv) for equiv in equivalences]...)
    allvec_t = Matrix{Float64}(permutedims(all_equiv_mat * lattvec))

    # Pre-compute CartesianIndex array for grid filling
    total_pts = sum(size(e, 1) for e in equivalences)
    all_indices = Vector{CartesianIndex{3}}(undef, total_pts)
    k = 1
    for equiv in equivalences
        for ir = 1:size(equiv, 1)
            r1 = equiv[ir, 1]
            r2 = equiv[ir, 2]
            r3 = equiv[ir, 3]
            all_indices[k] =
                CartesianIndex(mod(r1, d1) + 1, mod(r2, d2) + 1, mod(r3, d3) + 1)
            k += 1
        end
    end

    # Output arrays
    eband = zeros(nbands, npts)
    vvband = compute_velocity ? zeros(nbands, 3, 3, npts) : nothing

    # Call type-stable inner kernel
    _getBTPbands_kernel!(
        eband,
        vvband,
        coeffs,
        equivalences,
        allvec_t,
        all_indices,
        dims,
        plan,
    )

    return eband, vvband
end

#=
    getBTPbands_parallel(coeffs, equivalences, lattvec; compute_velocity=true)

Parallel version of getBTPbands using Julia threads.

Uses the same optimization strategy as getBTPbands:
1. Pre-computes Cartesian coordinates (lattvec * R)
2. Each thread gets its own work arrays
3. Fills all 4 grids in a single loop per band
=#
function getBTPbands_parallel(coeffs, equivalences, lattvec; compute_velocity = true)
    dims = determine_fft_dims(equivalences)
    nbands = size(coeffs, 1)
    npts = prod(dims)
    d1, d2, d3 = dims
    scale = prod(dims)

    # Pre-computation (amortized cost, shared across threads)
    plan = plan_ifft(zeros(ComplexF64, dims))
    all_equiv_mat = vcat([Matrix(equiv) for equiv in equivalences]...)
    allvec_t = Matrix{Float64}(permutedims(all_equiv_mat * lattvec))

    total_pts = sum(size(e, 1) for e in equivalences)
    all_indices = Vector{CartesianIndex{3}}(undef, total_pts)
    k = 1
    for equiv in equivalences
        for ir = 1:size(equiv, 1)
            r1 = equiv[ir, 1]
            r2 = equiv[ir, 2]
            r3 = equiv[ir, 3]
            all_indices[k] =
                CartesianIndex(mod(r1, d1) + 1, mod(r2, d2) + 1, mod(r3, d3) + 1)
            k += 1
        end
    end

    # Output arrays
    eband = zeros(nbands, npts)
    vvband = compute_velocity ? zeros(nbands, 3, 3, npts) : nothing

    # Parallel processing with thread-local work arrays
    Threads.@threads for iband = 1:nbands
        # Thread-local work arrays
        egrid = zeros(ComplexF64, dims)
        vgrid1 = zeros(ComplexF64, dims)
        vgrid2 = zeros(ComplexF64, dims)
        vgrid3 = zeros(ComplexF64, dims)

        flat_coeffs = _flatten_coeffs(coeffs[iband, :], equivalences)
        _fill_grids!(egrid, vgrid1, vgrid2, vgrid3, flat_coeffs, allvec_t, all_indices)

        eband[iband, :] = scale .* real.(vec(plan * egrid))

        if compute_velocity
            vb1 = scale .* real.(vec(plan * vgrid1))
            vb2 = scale .* real.(vec(plan * vgrid2))
            vb3 = scale .* real.(vec(plan * vgrid3))

            # Exploit symmetry: vvband[i,j] = vvband[j,i] (6 multiplications instead of 9)
            @inbounds for k = 1:npts
                # Diagonal components (3 multiplications)
                vvband[iband, 1, 1, k] = vb1[k] * vb1[k]
                vvband[iband, 2, 2, k] = vb2[k] * vb2[k]
                vvband[iband, 3, 3, k] = vb3[k] * vb3[k]

                # Off-diagonal components (3 multiplications + 6 assignments)
                v12 = vb1[k] * vb2[k]
                v13 = vb1[k] * vb3[k]
                v23 = vb2[k] * vb3[k]
                vvband[iband, 1, 2, k] = v12
                vvband[iband, 2, 1, k] = v12
                vvband[iband, 1, 3, k] = v13
                vvband[iband, 3, 1, k] = v13
                vvband[iband, 2, 3, k] = v23
                vvband[iband, 3, 2, k] = v23
            end
        end
    end

    return eband, vvband
end

#=
    getBTPbands_wannier(interp::WannierInterpolator, mesh::UniformMesh)

Wannier-path counterpart to `getBTPbands`. Generates a uniform mesh of
fractional k-points from `mesh.nk`, evaluates `interpolate_bands` and
`interpolate_velocities`, and constructs the velocity outer-product
tensor `vvband = v ⊗ v` so that the downstream `BTPDOS` / Onsager flow
in `solve_transport` is shared with the Fourier path.

Requires `mesh.nk` to be an explicit tuple; `UniformMesh(nothing)` errors
because the Wannier path has no equivalences-derived natural mesh.

Returns `(eband (n_wann, npts), vvband (n_wann, 3, 3, npts))` with energy
in Hartree and velocity outer products in (Hartree·Bohr)².
=#
# Energy gap (in Hartree) below which two interpolated bands are treated as a
# degenerate multiplet for the velocity-tensor construction. Equal to 1e-6 eV.
const WANNIER_DEGEN_THRESHOLD_HA = 1.0e-6 * EV_TO_HA

# Number of approach directions for the degenerate-manifold angular average.
const WANNIER_DEGEN_NDIRS = 200

# `n` quasi-uniform unit vectors on the sphere (Fibonacci lattice), used as the
# approach directions for the degenerate-manifold angular average.
function _fibonacci_directions(n::Int)
    dirs = Vector{NTuple{3,Float64}}(undef, n)
    ga = π * (3 - sqrt(5))
    for i = 0:(n-1)
        z = 1 - 2 * (i + 0.5) / n
        r = sqrt(max(0.0, 1 - z * z))
        θ = ga * i
        dirs[i+1] = (r * cos(θ), r * sin(θ), z)
    end
    return dirs
end

# Group band indices whose energies are degenerate within `thr` (Hartree).
# `energies` is the per-k band-energy vector; returns groups of indices into it.
function _degenerate_groups(energies::AbstractVector, thr::Real)
    n = length(energies)
    order = sortperm(energies)
    groups = Vector{Vector{Int}}()
    current = [order[1]]
    @inbounds for t = 2:n
        if energies[order[t]] - energies[order[t-1]] < thr
            push!(current, order[t])
        else
            push!(groups, current)
            current = [order[t]]
        end
    end
    push!(groups, current)
    return groups
end

# Write the rank-1 velocity outer product v ⊗ v into vvband for band `b` at `k`.
@inline function _fill_vv!(vvband, b, k, v1, v2, v3)
    vvband[b, 1, 1, k] = v1 * v1
    vvband[b, 2, 2, k] = v2 * v2
    vvband[b, 3, 3, k] = v3 * v3
    v12 = v1 * v2
    v13 = v1 * v3
    v23 = v2 * v3
    vvband[b, 1, 2, k] = v12
    vvband[b, 2, 1, k] = v12
    vvband[b, 1, 3, k] = v13
    vvband[b, 3, 1, k] = v13
    vvband[b, 2, 3, k] = v23
    vvband[b, 3, 2, k] = v23
    return nothing
end

# Relative eigenvalue-gap below which B(d) is treated as itself degenerate along an
# approach direction; such directions are skipped (their split basis would again be
# LAPACK-dependent). The cut is on eigenvalues, which are platform-stable, so the
# skip decision is reproducible.
const WANNIER_DIRECTION_GAP_REL = 1.0e-6

# Gauge-invariant velocity tensor for a degenerate multiplet via the directional
# (angular-average) limit. For each approach direction `d`, the multiplet splits
# along the eigenbasis of B(d) = Σ_a d_a A_a; the band-diagonal velocities in that
# "good" basis give the semiclassical v⊗v contribution, averaged over directions.
# Directions where B(d) is itself (near-)degenerate are skipped — there the split
# basis is not unique and would reintroduce the LAPACK gauge dependence — and the
# average is taken over the directions actually used. The averaged block total is
# distributed equally over the `L` members (they share one DOS bin). If every
# direction is degenerate (B(d) ∝ I for all d, pathological), fall back to the
# gauge-invariant trace tr(A_i A_j).
function _degenerate_vv!(vvband, A, grp, k, dirs)
    L = length(grp)
    A1 = A[grp, grp, 1, k]
    A2 = A[grp, grp, 2, k]
    A3 = A[grp, grp, 3, k]
    Ab = (A1, A2, A3)
    acc = zeros(3, 3)
    vt = Matrix{Float64}(undef, L, 3)
    used = 0
    for d in dirs
        B = d[1] .* A1 .+ d[2] .* A2 .+ d[3] .* A3
        F = eigen(Hermitian(B))
        vals = F.values
        spread = vals[end] - vals[1]
        mingap = spread
        for m = 2:L
            mingap = min(mingap, vals[m] - vals[m-1])
        end
        # Skip directions where B(d) is itself (near-)degenerate.
        spread > 0 && mingap < WANNIER_DIRECTION_GAP_REL * spread && continue
        used += 1
        W = F.vectors
        for a = 1:3
            M = W' * Ab[a] * W
            for m = 1:L
                vt[m, a] = real(M[m, m])
            end
        end
        for i = 1:3, j = 1:3
            s = 0.0
            for m = 1:L
                s += vt[m, i] * vt[m, j]
            end
            acc[i, j] += s
        end
    end
    if used == 0
        # Pathological: every direction is degenerate. Use the gauge-invariant
        # trace tr(A_i A_j) (also distributed over the block).
        for i = 1:3, j = 1:3
            acc[i, j] = real(tr(Ab[i] * Ab[j]))
        end
        used = 1
    end
    denom = used * L
    for i = 1:3, j = 1:3
        val = acc[i, j] / denom
        for b in grp
            vvband[b, i, j, k] = val
        end
    end
    return nothing
end

function getBTPbands_wannier(interp::WannierInterpolator, mesh::UniformMesh)
    isnothing(mesh.nk) && error(
        "WannierInterpolator requires an explicit UniformMesh size, got " *
        "UniformMesh(nothing). Pass UniformMesh((nk1, nk2, nk3)).",
    )
    n1, n2, n3 = mesh.nk
    npts = n1 * n2 * n3
    kpts = Matrix{Float64}(undef, npts, 3)
    idx = 1
    for k3 = 0:(n3-1), k2 = 0:(n2-1), k1 = 0:(n1-1)
        kpts[idx, 1] = k1 / n1
        kpts[idx, 2] = k2 / n2
        kpts[idx, 3] = k3 / n3
        idx += 1
    end

    eband = interpolate_bands(interp, kpts)
    # Full covariant velocity matrices (band gauge); the diagonal is the group
    # velocity, the intra-multiplet off-diagonal feeds the degenerate construction.
    A = velocity_matrices(interp, kpts)

    n_wann = size(eband, 1)
    vvband = zeros(n_wann, 3, 3, npts)
    dirs = _fibonacci_directions(WANNIER_DEGEN_NDIRS)
    for k = 1:npts
        groups = _degenerate_groups(view(eband, :, k), WANNIER_DEGEN_THRESHOLD_HA)
        for grp in groups
            if length(grp) == 1
                b = grp[1]
                _fill_vv!(
                    vvband,
                    b,
                    k,
                    real(A[b, b, 1, k]),
                    real(A[b, b, 2, k]),
                    real(A[b, b, 3, k]),
                )
            else
                _degenerate_vv!(vvband, A, grp, k, dirs)
            end
        end
    end

    return eband, vvband
end
