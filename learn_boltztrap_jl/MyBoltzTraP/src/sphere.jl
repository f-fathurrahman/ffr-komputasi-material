# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/sphere

using LinearAlgebra

# Note: SpinType (Unpolarized, Collinear, NonCollinear) is defined in spintypes.jl
# and included before this file. Used for magnetic moment type dispatch.

#=
    compute_bounds(lattvec::AbstractMatrix, radius::Real) -> Vector{Int}

Compute search bounds for lattice points within a sphere.

Given a radius r in reciprocal space, determine the maximum |n_i| for each
lattice direction such that all points with |R| ≤ r are included.

Uses the metric tensor G = L^T L to compute:
    bounds[k] = ceil(r × √(G⁻¹ₖₖ))

# Arguments
- `lattvec`: 3×3 matrix of lattice vectors (columns)
- `radius`: Sphere radius in reciprocal space units

# Returns
- Vector of 3 integers representing search bounds [n1_max, n2_max, n3_max]

# Example
```julia
lattvec = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]  # Simple cubic
bounds = compute_bounds(lattvec, 50.0)  # → [10, 10, 10]
```
=#
function compute_bounds(lattvec::AbstractMatrix, radius::Real)
    # Metric tensor G = L^T L
    metric = lattvec' * lattvec

    # Inverse metric G^{-1}
    invmetric = inv(metric)

    # bounds[k] = ceil(r × √(G⁻¹ₖₖ))
    bounds = Vector{Int}(undef, 3)
    for k = 1:3
        bounds[k] = ceil(Int, radius * sqrt(invmetric[k, k]))
    end

    return bounds
end

#=
    compute_radius(lattvec, nrotations, nkpt_target) -> Float64

Estimate the sphere radius needed to obtain approximately `nkpt_target`
equivalence classes.

Port of BoltzTraP2/sphere/__init__.py compute_radius function.

The formula is based on the volume of a sphere in reciprocal space:
    npoints = nkpt_target × nrotations
    radius = (3/(4π) × npoints × vol)^(1/3)

# Arguments
- `lattvec`: 3×3 matrix of lattice vectors (columns)
- `nrotations`: Number of unique rotation operations
- `nkpt_target`: Target number of equivalence classes (k-points)

# Returns
- Estimated sphere radius in reciprocal space units

# Example
```julia
lattvec = 5.43 * [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]  # Si FCC
nrot = 48
radius = compute_radius(lattvec, nrot, 5000)  # ≈ 127.8
```
=#
function compute_radius(lattvec::AbstractMatrix, nrotations::Integer, nkpt_target::Integer)
    vol = abs(det(lattvec))
    npoints = nkpt_target * nrotations
    return (3.0 / (4π) * npoints * vol)^(1 / 3)
end

#=
    find_permutation(lattvec, positions, types, rotation, translation, symprec)

Find the atomic permutation induced by a symmetry operation.

Returns a vector `perm` where `perm[i]` is the index of the atom that
atom `i` maps to under the symmetry operation.

# Arguments
- `lattvec`: 3×3 lattice vectors (columns)
- `positions`: Vector of fractional positions (each is 3-element vector)
- `types`: Vector of atom type indices
- `rotation`: 3×3 rotation matrix (in fractional coordinates)
- `translation`: 3-element translation vector (in fractional coordinates)
- `symprec`: Symmetry precision tolerance
=#
function find_permutation(
    lattvec::AbstractMatrix,
    positions::AbstractVector,
    types::AbstractVector{<:Integer},
    rotation::AbstractMatrix,
    translation::AbstractVector,
    symprec::Real,
)
    natoms = length(positions)
    perm = Vector{Int}(undef, natoms)

    for i = 1:natoms
        # Apply symmetry operation: R * pos + t
        pos_i = positions[i]
        new_pos = rotation * pos_i + translation

        # Wrap to [0, 1) range
        new_pos = mod.(new_pos, 1.0)

        # Find matching atom
        found = false
        for j = 1:natoms
            if types[i] != types[j]
                continue
            end

            # Compute difference in fractional coordinates
            diff = new_pos - positions[j]
            # Account for periodic boundary
            diff = mod.(diff .+ 0.5, 1.0) .- 0.5

            # Convert to Cartesian and check distance
            cart_diff = lattvec * diff
            if norm(cart_diff) < symprec
                perm[i] = j
                found = true
                break
            end
        end

        if !found
            error("Could not find permutation for atom $i")
        end
    end

    return perm
end

#=
    determine_compatibility(rotation, perm, magmom, spintype, symprec)

Check if a symmetry operation is compatible with the magnetic configuration.

Returns a named tuple `(forward=Bool, backward=Bool)` indicating if the
forward operation and time-reversal operation are compatible.

# Arguments
- `rotation`: 3×3 Cartesian rotation matrix
- `perm`: Atomic permutation vector
- `magmom`: Magnetic moments (`nothing`, `Vector{Float64}`, or `Vector` of 3-vectors)
- `spintype`: SpinType instance (Unpolarized(), Collinear(), or NonCollinear())
- `symprec`: Symmetry precision
=#
function determine_compatibility(
    rotation::AbstractMatrix,
    perm::AbstractVector{<:Integer},
    magmom,
    spintype::SpinType,
    symprec::Real,
)
    if spintype isa Unpolarized
        # No magnetic constraints - both forward and backward are compatible
        return (forward = true, backward = true)
    elseif spintype isa Collinear
        # Collinear: check if moments match after permutation
        # Allow for spin flip (time reversal)
        natoms = length(perm)
        forward_ok = true
        backward_ok = true

        for i = 1:natoms
            j = perm[i]
            m_i = magmom[i]
            m_j = magmom[j]

            # Forward: moments should match (|m_i - m_j| < symprec)
            if abs(m_i - m_j) >= symprec
                forward_ok = false
            end

            # Backward (time reversal): moments can flip (|m_i + m_j| < symprec)
            if abs(m_i + m_j) >= symprec
                backward_ok = false
            end
        end

        return (forward = forward_ok, backward = backward_ok)
    else  # NonCollinear
        # Non-collinear: apply rotation to magnetic moment vectors
        natoms = length(perm)
        forward_ok = true
        backward_ok = true

        for i = 1:natoms
            j = perm[i]
            m_i = magmom[i]
            m_j = magmom[j]

            # Forward: R * m_i should equal m_j
            rotated = rotation * m_i
            if norm(rotated - m_j) >= symprec
                forward_ok = false
            end

            # Backward (time reversal): R * m_i should equal -m_j
            if norm(rotated + m_j) >= symprec
                backward_ok = false
            end
        end

        return (forward = forward_ok, backward = backward_ok)
    end
end

#=
    get_spintype(magmom) -> SpinType

Determine the spin type from the magnetic moment input.

- `nothing` → Unpolarized()
- `Vector{<:Real}` → Collinear()
- `Vector{<:AbstractVector}` (3D vectors) → NonCollinear()
=#
function get_spintype(magmom)
    if isnothing(magmom)
        return Unpolarized()
    elseif eltype(magmom) <: Real
        return Collinear()
    else
        return NonCollinear()
    end
end

# Backward compatibility alias
const get_magmom_type = get_spintype

#=
    calc_nrotations(lattvec, positions, types, magmom; symprec=1e-5)

Compute the number of unique rotations taking time reversal and magnetic
constraints into account.

# Arguments
- `lattvec`: 3×3 lattice vectors (columns)
- `positions`: Vector of fractional positions (N×3 or Vector of 3-vectors)
- `types`: Vector of atom type indices
- `magmom`: Magnetic moments (`nothing` for unpolarized)
- `symprec`: Symmetry precision

# Returns
- Number of unique rotation operations
=#
function calc_nrotations(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Convert positions to vector of SVectors for Spglib
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:size(positions, 1)]
    return calc_nrotations(lattvec, pos_vec, types, magmom; symprec)
end

function calc_nrotations(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Determine magnetic type
    mtype = get_magmom_type(magmom)

    # Create Spglib cell
    cell = Spglib.Cell(lattvec, positions, types)

    # Get symmetry dataset
    dataset = Spglib.get_dataset(cell, symprec)

    # Get rotations and translations
    rotations = dataset.rotations
    translations = dataset.translations

    # Compute unique rotations with time-reversal
    unique_rots = Set{Matrix{Int}}()

    # Prepare for Cartesian conversion
    # Correct formula: R_cart = L * R_frac * L^{-1}
    L = lattvec
    L_inv = inv(L)

    for (iop, rot) in enumerate(rotations)
        trans = translations[iop]

        # Get atomic permutation
        perm = find_permutation(lattvec, positions, types, rot, trans, symprec)

        # Convert rotation to Cartesian coordinates
        rot_cart = L * Float64.(rot) * L_inv

        # Check magnetic compatibility
        compat = determine_compatibility(rot_cart, perm, magmom, mtype, symprec)

        # Store rotation directly (NOT transposed)
        # The C++ code stores wrapper.transpose() where wrapper is already R^T
        # due to Eigen's column-major interpretation of row-major spglib data,
        # so C++ stores (R^T)^T = R. We should store R directly.
        rot_int = Matrix{Int}(rot)

        if compat.forward
            push!(unique_rots, rot_int)
        end
        if compat.backward
            # Time reversal: -R
            push!(unique_rots, -rot_int)
        end
    end

    return length(unique_rots)
end

#=
    get_spacegroup_info(lattvec, positions, types; symprec=1e-5)

Get space group information from crystal structure.

# Returns
- `spacegroup_number`: International space group number (1-230)
- `international_symbol`: International symbol (e.g., "Fd-3m")
=#
function get_spacegroup_info(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:size(positions, 1)]
    return get_spacegroup_info(lattvec, pos_vec, types; symprec)
end

function get_spacegroup_info(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    cell = Spglib.Cell(lattvec, positions, types)
    dataset = Spglib.get_dataset(cell, symprec)
    return (
        spacegroup_number = dataset.spacegroup_number,
        international_symbol = String(dataset.international_symbol),
    )
end

#=
    get_unique_rotations(lattvec, positions, types, magmom; symprec=1e-5)

Get unique rotation matrices taking time reversal and magnetic constraints
into account. Returns both the integer rotations and Cartesian rotations.

Note: The C++ BoltzTraP2 code stores `wrapper.transpose()` where `wrapper` is
an Eigen::Map of spglib's row-major rotation data. Since Eigen uses column-major
by default, `wrapper` is actually R^T, so `wrapper.transpose()` = R.
We store R directly (not R^T) to match the C++ behavior.

# Returns
- `rotations`: Vector of 3×3 integer rotation matrices (for lattice point operations)
- `crotations`: Vector of 3×3 Cartesian rotation matrices
=#
function get_unique_rotations(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Determine magnetic type
    mtype = get_magmom_type(magmom)

    # Create Spglib cell
    cell = Spglib.Cell(lattvec, positions, types)

    # Get symmetry dataset
    dataset = Spglib.get_dataset(cell, symprec)

    # Get rotations and translations
    rotations = dataset.rotations
    translations = dataset.translations

    # Compute unique rotations with time-reversal
    unique_rots = Set{Matrix{Int}}()
    rot_to_crot = Dict{Matrix{Int},Matrix{Float64}}()

    # Prepare for Cartesian conversion
    # Correct formula: R_cart = L * R_frac * L^{-1}
    L = lattvec
    L_inv = inv(L)

    for (iop, rot) in enumerate(rotations)
        trans = translations[iop]

        # Get atomic permutation
        perm = find_permutation(lattvec, positions, types, rot, trans, symprec)

        # Convert rotation to Cartesian coordinates: R_cart = L * R * L^{-1}
        rot_cart = L * Float64.(rot) * L_inv

        # Check magnetic compatibility
        compat = determine_compatibility(rot_cart, perm, magmom, mtype, symprec)

        # Store rotation directly (NOT transposed)
        # The C++ code stores wrapper.transpose() where wrapper is already R^T
        # due to Eigen's column-major interpretation of row-major spglib data,
        # so C++ stores (R^T)^T = R. We should store R directly.
        rot_int = Matrix{Int}(rot)

        if compat.forward
            if rot_int ∉ unique_rots
                push!(unique_rots, rot_int)
                # Cartesian rotation: R_cart = L * R * L^{-1}
                rot_to_crot[rot_int] = rot_cart
            end
        end
        if compat.backward
            neg_rot = -rot_int
            if neg_rot ∉ unique_rots
                push!(unique_rots, neg_rot)
                # Cartesian rotation for backward (time reversal): -R_cart
                rot_to_crot[neg_rot] = -rot_cart
            end
        end
    end

    # Convert to vectors
    rots_vec = collect(unique_rots)
    crots_vec = [rot_to_crot[r] for r in rots_vec]

    return rots_vec, crots_vec
end

function get_unique_rotations(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:size(positions, 1)]
    return get_unique_rotations(lattvec, pos_vec, types, magmom; symprec)
end

#=
    calc_tensor_basis(lattvec, positions, types, magmom; symprec=1e-5)

Compute a basis for 3×3 tensors compatible with the crystal symmetry.

Returns a set of independent 3×3 tensors that span all possible 3×3
tensors compatible with the symmetries of the system. Any linear response
tensor (like conductivity) must be a linear combination of these basis tensors.

# Arguments
- `lattvec`: 3×3 lattice vectors (columns)
- `positions`: Atomic positions (N×3 matrix or Vector of 3-vectors)
- `types`: Vector of atom type indices
- `magmom`: Magnetic moments (`nothing` for unpolarized)
- `symprec`: Symmetry precision

# Returns
- Array of shape (nbasis, 3, 3) containing the basis tensors
=#
function calc_tensor_basis(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:size(positions, 1)]
    return calc_tensor_basis(lattvec, pos_vec, types, magmom; symprec)
end

function calc_tensor_basis(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom;
    symprec::Real = 1e-5,
)
    # Get unique Cartesian rotations
    _, crotations = get_unique_rotations(lattvec, positions, types, magmom; symprec)

    nops = length(crotations)

    # Build constraint matrix
    # For each rotation R, the constraint is: R^T T R = T (tensor invariance)
    # In component form: R_βα R_γδ T_βγ = T_αδ
    # Rearranged: (R_βα R_γδ - δ_βα δ_γδ) T_βγ = 0
    #
    # Total constraints: 9 per rotation + 3 for symmetry
    constraints = zeros(9 * nops + 3, 9)

    I3 = Matrix{Float64}(I, 3, 3)

    for (iop, rot) in enumerate(crotations)
        offset = 9 * (iop - 1)

        for α = 1:3
            for δ = 1:3
                row = offset + 3 * (α - 1) + δ
                for β = 1:3
                    for γ = 1:3
                        col = 3 * (β - 1) + γ
                        constraints[row, col] = rot[β, α] * rot[γ, δ] - I3[β, α] * I3[γ, δ]
                    end
                end
            end
        end
    end

    # Add symmetry constraints: T_ij = T_ji
    # T[0,1] = T[1,0]: col 2 (0*3+1+1) - col 4 (1*3+0+1)
    # T[0,2] = T[2,0]: col 3 (0*3+2+1) - col 7 (2*3+0+1)
    # T[1,2] = T[2,1]: col 6 (1*3+2+1) - col 8 (2*3+1+1)
    sym_offset = 9 * nops
    constraints[sym_offset+1, 2] = 1.0   # T[1,2] in 1-indexed = T[0,1] in 0-indexed
    constraints[sym_offset+1, 4] = -1.0  # T[2,1] in 1-indexed = T[1,0] in 0-indexed
    constraints[sym_offset+2, 3] = 1.0   # T[1,3] = T[0,2]
    constraints[sym_offset+2, 7] = -1.0  # T[3,1] = T[2,0]
    constraints[sym_offset+3, 6] = 1.0   # T[2,3] = T[1,2]
    constraints[sym_offset+3, 8] = -1.0  # T[3,2] = T[2,1]

    # Compute null space using SVD
    U, S, V = svd(constraints)

    # Find null space vectors (singular values ~ 0)
    tol = 1e-10 * max(size(constraints)...) * maximum(S)
    null_indices = findall(s -> s < tol, S)

    # The null space is spanned by the columns of V corresponding to zero singular values
    # For SVD of A, the null space of A is spanned by columns of V where σ ≈ 0
    nbasis = length(null_indices)

    if nbasis == 0
        # Check if there are very small singular values
        # The null space dimension should be 9 - rank(constraints)
        rank_approx = count(s -> s >= tol, S)
        nbasis = 9 - rank_approx
        if nbasis > 0
            null_indices = (9-nbasis+1):9
        end
    end

    # Extract basis tensors from null space
    if nbasis == 0
        return zeros(0, 3, 3)
    end

    basis = zeros(nbasis, 3, 3)
    for (i, idx) in enumerate(null_indices)
        vec = V[:, idx]
        # Reshape to 3×3 (row-major as in C++)
        for r = 1:3
            for c = 1:3
                basis[i, r, c] = vec[3*(r-1)+c]
            end
        end
    end

    return basis
end
