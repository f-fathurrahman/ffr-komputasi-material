# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/sphere

#=
    lattice_points_in_sphere(lattvec, radius, bounds)

Enumerate all lattice points within a sphere of given radius.

# Arguments
- `lattvec`: 3×3 lattice vectors (columns)
- `radius`: Sphere radius in reciprocal space
- `bounds`: (n1, n2, n3) search bounds for each direction

# Returns
- `points`: Vector of SVector{3,Int} lattice points, sorted by squared norm
- `sqnorms`: Vector of squared norms corresponding to each point
=#
function lattice_points_in_sphere(lattvec, radius, bounds)
    metric = lattvec' * lattvec  # G = L^T L
    r2 = radius^2
    points = Vector{SVector{3,Int}}()
    sqnorms = Vector{Float64}()

    for I in CartesianIndices((
        (-bounds[1]):bounds[1],
        (-bounds[2]):bounds[2],
        (-bounds[3]):bounds[3],
    ))
        n = SVector{3,Int}(Tuple(I))
        norm_sq = n' * metric * n
        if norm_sq < r2  # Strict inequality, matching C++ code
            push!(points, n)
            push!(sqnorms, norm_sq)
        end
    end

    # Sort by squared norm
    perm = sortperm(sqnorms)
    return points[perm], sqnorms[perm]
end

#=
    compute_equivalence_classes(points, sqnorms, rotations, symprec=1e-5)

Compute equivalence classes for lattice points under symmetry operations.

Uses the same algorithm as C++ BoltzTraP2: for each point, only compare with
earlier points that have the same squared norm (within tolerance). This is
much faster than checking all points since equivalent points must have the
same norm.

# Arguments
- `points`: Vector of lattice points (sorted by squared norm)
- `sqnorms`: Vector of squared norms corresponding to each point
- `rotations`: Vector of rotation matrices
- `symprec`: Tolerance for norm comparison

# Returns
- `mapping`: Vector where mapping[i] gives the representative index for point i
=#
function compute_equivalence_classes(points, sqnorms, rotations, symprec = 1e-5)
    n = length(points)

    # Build a lookup table for points
    point_to_idx = Dict{SVector{3,Int},Int}()
    for (i, p) in enumerate(points)
        point_to_idx[p] = i
    end

    # Mapping: -1 means not yet assigned
    mapping = fill(-1, n)

    for i = 1:n
        point = points[i]
        sqnorm = sqnorms[i]

        # Points can only be equivalent if they have the same norm
        # Binary search for range of candidates with similar norm
        lo = (1.0 - symprec) * sqnorm
        hi = (1.0 + symprec) * sqnorm
        lbound = searchsortedfirst(sqnorms, lo)
        ubound = min(searchsortedlast(sqnorms, hi), i - 1)  # Only earlier points

        # Check each rotation
        for R in rotations
            image = SVector{3,Int}(R * point)
            j = get(point_to_idx, image, 0)
            if j > 0 && lbound <= j <= ubound
                # Found a match - find its representative
                k = j
                while mapping[k] != k
                    k = mapping[k]
                end
                mapping[i] = k
                break
            end
        end

        if mapping[i] != -1
            # Already found a class, skip remaining rotations
            continue
        end

        # No equivalence found - start a new class
        if mapping[i] == -1
            mapping[i] = i
        end
    end

    return mapping
end

#=
    group_by_equivalence(points, mapping)

Group points into their equivalence classes.

# Returns
- Vector of matrices, where each matrix contains all equivalent points (as rows)
=#
function group_by_equivalence(points, mapping)
    classes = Dict{Int,Vector{SVector{3,Int}}}()
    for (i, rep) in enumerate(mapping)
        if !haskey(classes, rep)
            classes[rep] = SVector{3,Int}[]
        end
        push!(classes[rep], points[i])
    end

    # Convert to vector of matrices, sorted by representative index
    reps = sort(collect(keys(classes)))
    return [hcat(classes[r]...)' for r in reps]
end

#=
    calc_sphere_quotient_set(lattvec, rotations, radius, bounds; symprec=1e-5)

Compute equivalence classes for lattice points in a sphere using given rotations.

# Returns
- `equivalences`: Vector of matrices, each containing equivalent lattice points (as rows)
=#
function calc_sphere_quotient_set(lattvec, rotations, radius, bounds; symprec = 1e-5)
    points, sqnorms = lattice_points_in_sphere(lattvec, radius, bounds)
    mapping = compute_equivalence_classes(points, sqnorms, rotations, symprec)
    return group_by_equivalence(points, mapping)
end

#=
    calc_sphere_quotient_set(lattvec, positions, types, magmom, radius, bounds; symprec=1e-5)

Compute equivalence classes for lattice points in a sphere.

This is the main entry point, equivalent to Python's calc_sphere_quotient_set.
It automatically computes the symmetry rotations from the crystal structure.

# Arguments
- `lattvec`: 3×3 lattice vectors (columns)
- `positions`: Atomic positions (N×3 matrix or Vector of 3-vectors)
- `types`: Vector of atom type indices
- `magmom`: Magnetic moments (`nothing` for unpolarized)
- `radius`: Sphere radius in reciprocal space
- `bounds`: Search bounds [n1_max, n2_max, n3_max]
- `symprec`: Symmetry precision

# Returns
- `equivalences`: Vector of vectors, where each inner vector contains
  equivalent lattice points as SVector{3,Int}
=#
function calc_sphere_quotient_set(
    lattvec::AbstractMatrix,
    positions::AbstractMatrix,
    types::AbstractVector{<:Integer},
    magmom,
    radius::Real,
    bounds::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    pos_vec = [SVector{3,Float64}(positions[i, :]) for i = 1:size(positions, 1)]
    return calc_sphere_quotient_set(
        lattvec,
        pos_vec,
        types,
        magmom,
        radius,
        bounds;
        symprec,
    )
end

function calc_sphere_quotient_set(
    lattvec::AbstractMatrix,
    positions::AbstractVector{<:AbstractVector},
    types::AbstractVector{<:Integer},
    magmom,
    radius::Real,
    bounds::AbstractVector{<:Integer};
    symprec::Real = 1e-5,
)
    # Get unique rotations (integer matrices for reciprocal space)
    rotations, _ = get_unique_rotations(lattvec, positions, types, magmom; symprec)

    # Compute equivalence classes
    points, sqnorms = lattice_points_in_sphere(lattvec, radius, bounds)
    mapping = compute_equivalence_classes(points, sqnorms, rotations, symprec)

    # Group by equivalence - return as Vector of Vectors (easier to compare)
    classes = Dict{Int,Vector{SVector{3,Int}}}()
    for (i, rep) in enumerate(mapping)
        if !haskey(classes, rep)
            classes[rep] = SVector{3,Int}[]
        end
        push!(classes[rep], points[i])
    end

    # Sort by representative index to get consistent ordering
    reps = sort(collect(keys(classes)))
    return [classes[r] for r in reps]
end
