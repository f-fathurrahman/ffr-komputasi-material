# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl
"""
Type definitions for BoltzTraP.jl.
"""

"""
    DFTData{NSpin}

Container for DFT calculation results with spin parameter.

# Type Parameter

  - `NSpin::Int`: Number of spin channels

      + `1`: Non-spin-polarized (non-magnetic)
      + `2`: Spin-polarized (collinear magnetic, ISPIN=2)

# Fields

  - `lattice::Matrix{Float64}`: Lattice vectors (3×3) in Bohr, columns are vectors
  - `positions::Matrix{Float64}`: Atomic positions (3×natom) in fractional coordinates
  - `species::Vector{String}`: Element symbols
  - `kpoints::Matrix{Float64}`: K-points (3×nkpts) in fractional coordinates
  - `weights::Vector{Float64}`: K-point weights (nkpts,)
  - `ebands::Array{Float64,3}`: Band energies (nbands, nkpts, nspin) in Hartree
  - `occupations::Array{Float64,3}`: Occupations (nbands, nkpts, nspin)
  - `fermi::Float64`: Fermi energy in Hartree
  - `nelect::Float64`: Number of electrons
  - `magmom::Union{Nothing,Vector{Float64},Matrix{Float64}}`: Magnetic moments (optional)

# Notes

  - Non-collinear (spinor) calculations are not supported in v0.1
  - When `NSpin=2`, bands from both spin channels are stored separately
"""
struct DFTData{NSpin}
    lattice::Matrix{Float64}
    positions::Matrix{Float64}
    species::Vector{String}
    kpoints::Matrix{Float64}
    weights::Vector{Float64}
    ebands::Array{Float64,3}
    occupations::Array{Float64,3}
    fermi::Float64
    nelect::Float64
    magmom::Union{Nothing,Vector{Float64},Matrix{Float64}}
end

# Type aliases for clarity
const NonMagneticData = DFTData{1}
const SpinPolarizedData = DFTData{2}

#=
    DFTData(; lattice, positions, species, kpoints, weights,
              ebands, occupations, fermi, nelect, magmom=nothing) -> DFTData{NSpin}

Construct DFTData with NSpin inferred from ebands shape.
=#
function DFTData(;
    lattice::AbstractMatrix,
    positions::AbstractMatrix,
    species::AbstractVector{<:AbstractString},
    kpoints::AbstractMatrix,
    weights::AbstractVector,
    ebands::AbstractArray{<:Real,3},
    occupations::AbstractArray{<:Real,3},
    fermi::Real,
    nelect::Real,
    magmom::Union{Nothing,AbstractVector,AbstractMatrix} = nothing,
)
    nspin = size(ebands, 3)
    DFTData{nspin}(
        Matrix{Float64}(lattice),
        Matrix{Float64}(positions),
        Vector{String}(species),
        Matrix{Float64}(kpoints),
        Vector{Float64}(weights),
        Array{Float64,3}(ebands),
        Array{Float64,3}(occupations),
        Float64(fermi),
        Float64(nelect),
        isnothing(magmom) ? nothing :
        (magmom isa AbstractMatrix ? Matrix{Float64}(magmom) : Vector{Float64}(magmom)),
    )
end

"""
    nspin(data::DFTData{N}) -> Int

Return the number of spin channels.
"""
nspin(::DFTData{N}) where {N} = N

"""
    is_magnetic(data::DFTData) -> Bool

Return `true` if data is from spin-polarized calculation.
"""
is_magnetic(data::DFTData) = nspin(data) > 1

# ============================================================================
# BandData{ST} - Parametric band data container for v0.3 magnetic support
# ============================================================================

"""
    BandData{ST<:SpinType}

Container for band structure data with SpinType parameter.

# Type Parameter

  - `ST<:SpinType`: Spin polarization type
      + `Unpolarized`: Non-spin-polarized (dosweight=2.0, nspinor=1)
      + `Collinear`: Collinear spin-polarized (dosweight=1.0, nspinor=2)
      + `NonCollinear`: Non-collinear spin-polarized (dosweight=1.0, nspinor=1)

# Fields

  - `lattice::Matrix{Float64}`: Lattice vectors (3×3) in Bohr, columns are vectors
  - `positions::Matrix{Float64}`: Atomic positions (3×natom) in fractional coordinates
  - `species::Vector{Symbol}`: Element symbols
  - `kpoints::Matrix{Float64}`: K-points (3×nkpts) in fractional coordinates
  - `ebands::Array{Float64}`: Band energies in Hartree
      + Unpolarized/NonCollinear: (nbands, nkpts)
      + Collinear: (nbands, nkpts, 2)
  - `fermi::Float64`: Fermi energy in Hartree
  - `nelect::Float64`: Number of electrons
  - `magmom`: Magnetic moments
      + `nothing` for Unpolarized
      + `Vector{Float64}` for Collinear (scalar per atom)
      + `Vector{Vector{Float64}}` for NonCollinear (3D vector per atom)

# Example

```julia
# Unpolarized
data = BandData(;
    lattice = L,
    positions = P,
    species = S,
    kpoints = K,
    ebands = E,
    fermi = 0.2,
    nelect = 8.0,
)

# Collinear
data = BandData(;
    lattice = L,
    positions = P,
    species = S,
    kpoints = K,
    ebands = E,
    fermi = 0.2,
    nelect = 8.0,
    magmom = [1.0, -1.0],
)
```

Construct BandData with SpinType inferred from magmom.
See also: [`SpinType`](@ref), [`Unpolarized`](@ref), [`Collinear`](@ref), [`NonCollinear`](@ref)
"""
struct BandData{ST<:SpinType}
    lattice::Matrix{Float64}
    positions::Matrix{Float64}
    species::Vector{Symbol}
    kpoints::Matrix{Float64}
    ebands::Array{Float64}
    fermi::Float64
    nelect::Float64
    magmom::Union{Nothing,Vector{Float64},Vector{Vector{Float64}}}
end

#=
    BandData(; lattice, positions, species, kpoints, ebands, fermi, nelect,
               magmom=nothing) -> BandData{ST}

Construct BandData with SpinType inferred from magmom.

SpinType inference:
- `magmom === nothing` → `Unpolarized`
- `magmom isa Vector{<:Real}` → `Collinear`
- `magmom isa Vector{<:AbstractVector}` → `NonCollinear`
=#
function BandData(;
    lattice::AbstractMatrix,
    positions::AbstractMatrix,
    species::AbstractVector{<:Union{Symbol,AbstractString}},
    kpoints::AbstractMatrix,
    ebands::AbstractArray{<:Real},
    fermi::Real,
    nelect::Real,
    magmom = nothing,
)
    # Convert to standard types
    lat = Matrix{Float64}(lattice)
    pos = Matrix{Float64}(positions)
    spec = [s isa Symbol ? s : Symbol(s) for s in species]
    kpts = Matrix{Float64}(kpoints)
    bands = Array{Float64}(ebands)
    f = Float64(fermi)
    ne = Float64(nelect)

    # Determine SpinType and convert magmom
    if isnothing(magmom)
        return BandData{Unpolarized}(lat, pos, spec, kpts, bands, f, ne, nothing)
    elseif eltype(magmom) <: Real
        # Collinear: scalar magmom per atom
        mm = Vector{Float64}(magmom)
        return BandData{Collinear}(lat, pos, spec, kpts, bands, f, ne, mm)
    else
        # NonCollinear: 3D vector magmom per atom
        mm = [Vector{Float64}(m) for m in magmom]
        return BandData{NonCollinear}(lat, pos, spec, kpts, bands, f, ne, mm)
    end
end

# Error dispatch for invalid argument types
function BandData(x::Any)
    throw(
        ArgumentError(
            "BandData requires keyword arguments (lattice, positions, species, kpoints, ebands, fermi, nelect), got $(typeof(x))",
        ),
    )
end

"""
    get_dosweight(data::BandData{ST}) -> Float64

Return the DOS weight (spin degeneracy factor) for the band data.
"""
get_dosweight(::BandData{ST}) where {ST<:SpinType} = dosweight(ST)

"""
    get_nspinor(data::BandData{ST}) -> Int

Return the number of spin channels for the band data.
"""
get_nspinor(::BandData{ST}) where {ST<:SpinType} = nspinor(ST)
