# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/sphere

#=
    get_rotations(cell; magmoms=nothing, symprec=1e-5)

Get symmetry rotation matrices for a crystal structure.

# Arguments
- `cell`: Spglib cell (lattice, positions, types)
- `magmoms`: Magnetic moments (`nothing`, `Vector{Float64}`, or `Vector{SVector{3,Float64}}`)
- `symprec`: Symmetry precision

# Returns
- `rotations`: Vector of 3×3 integer rotation matrices
=#
function get_rotations(cell; magmoms = nothing, symprec = 1e-5)
    # TODO: Implement using Spglib.jl
    # For non-magnetic: get_dataset(cell, symprec).rotations
    # For collinear: get_symmetry_with_collinear_spin
    # For noncollinear: custom implementation needed
    error("Not implemented yet")
end

#=
    to_reciprocal_rotations(rotations, lattvec)

Convert real-space rotations to reciprocal-space rotations.

In reciprocal space: k' = (R^T)^{-1} k
Using metric tensor: k' = G^{-1} R^T G k
where G = L^T L is the metric tensor.
=#
function to_reciprocal_rotations(rotations, lattvec)
    metric = lattvec' * lattvec
    metric_inv = inv(metric)
    return [metric_inv * R' * metric for R in rotations]
end

# Note: determine_compatibility is fully implemented in sphere.jl
# This stub is kept for backward compatibility but delegates to sphere.jl implementation
