# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl - Port of BoltzTraP2 GENE/Generic loader

#=
GENE/Generic file format loader.

The GENE format is a generic format used by BoltzTraP2 for storing band structure
data. Files use .structure and .energy extensions.

File formats:
- .structure: Crystal structure (lattice vectors + atomic positions)
- .energy: Band energies and optional velocities
=#

using LinearAlgebra

#=
    _ffloat_gene(s)

Parse a floating-point number from a string.
Handles Fortran D-notation if present.
=#
function _ffloat_gene(s::AbstractString)
    return parse(Float64, replace(lowercase(strip(s)), "d" => "e"))
end

#=
    _read_gene_struct(filename)

Read GENE .structure file.

Format:
- Line 1: Title (ignored)
- Lines 2-4: Lattice vectors (3×3, in Bohr, row-major)
- Line 5: Number of atoms
- Lines 6+: Element symbol + Cartesian coordinates (Bohr)

Returns NamedTuple with:
- lattice: 3×3 matrix (columns are lattice vectors, Bohr)
- positions: 3×natom matrix (fractional coordinates)
- species: Vector of element symbols
=#
function _read_gene_struct(filename::String)
    lines = readlines(filename)

    # Line 1: Title (ignored)
    # Lines 2-4: Lattice vectors (row-major in file → column-major in Julia)
    lattice = zeros(3, 3)
    for i = 1:3
        parts = split(lines[i+1])
        for j = 1:3
            lattice[j, i] = _ffloat_gene(parts[j])
        end
    end

    # Line 5: Number of atoms
    natom = parse(Int, split(lines[5])[1])

    # Lines 6+: Element + Cartesian coordinates (Bohr)
    species = String[]
    cart_positions = zeros(3, natom)
    for i = 1:natom
        parts = split(lines[5+i])
        push!(species, parts[1])
        for j = 1:3
            cart_positions[j, i] = _ffloat_gene(parts[j+1])
        end
    end

    # Convert Cartesian to fractional coordinates
    # frac = inv(lattice) * cart
    positions = lattice \ cart_positions

    return (lattice = lattice, positions = positions, species = species)
end

#=
    _read_gene_energy(filename)

Read GENE .energy file.

Format:
- Line 1: Title (ignored)
- Line 2: nk nspin efermi(Ry)
- For each spin channel:
  - For each k-point:
    - Header: kx ky kz nband
    - Band lines: energy [vx vy vz] (Rydberg units)

For spin-polarized calculations (nspin=2):
- K-points are duplicated for each spin channel
- Bands are concatenated (spin channels combined)

Returns NamedTuple with:
- fermi: Fermi energy (Hartree)
- dosweight: DOS weight (2.0 for non-spin, 1.0 for spin)
- kpoints: 3×nk matrix (fractional)
- ebands: nbands×nk matrix (Hartree, spin channels concatenated)
- mommat: Optional momentum matrix or `nothing`
=#
function _read_gene_energy(filename::String)
    lines = readlines(filename)

    # Line 1: Title (ignored)
    # Line 2: nk nspin efermi(Ry)
    header = split(lines[2])
    nk = parse(Int, header[1])
    nspin = parse(Int, header[2])
    efermi_ry = _ffloat_gene(header[3])

    # Read all k-point blocks
    linenumber = 2  # Start after header (1-indexed)
    minband = typemax(Int)
    ebands_list = Vector{Float64}[]  # Each element is a vector of band energies for one k-point
    mommat_list = Vector{Vector{Float64}}[]  # Each element is a vector of velocity vectors for one k-point
    kpoints_all = Vector{Float64}[]

    for ispin = 1:nspin
        for ik = 1:nk
            linenumber += 1
            kline = split(lines[linenumber])
            kx = _ffloat_gene(kline[1])
            ky = _ffloat_gene(kline[2])
            kz = _ffloat_gene(kline[3])
            nband = parse(Int, kline[4])

            if nband < minband
                minband = nband
            end

            push!(kpoints_all, [kx, ky, kz])

            eband = Float64[]
            vband = Vector{Float64}[]
            for ib = 1:nband
                linenumber += 1
                fields = split(lines[linenumber])
                e = _ffloat_gene(fields[1])
                push!(eband, e)

                if length(fields) == 4
                    v = [
                        _ffloat_gene(fields[2]),
                        _ffloat_gene(fields[3]),
                        _ffloat_gene(fields[4]),
                    ]
                    push!(vband, v)
                else
                    push!(vband, Float64[])
                end
            end
            push!(ebands_list, eband)
            push!(mommat_list, vband)
        end
    end

    # K-points: first nk only (redundant for spin channels)
    kpoints = zeros(3, nk)
    for ik = 1:nk
        kpoints[:, ik] = kpoints_all[ik]
    end

    # Bands: concatenate spin channels
    # Shape: (nspin * minband, nk)
    nbands_total = nspin * minband
    ebands = zeros(nbands_total, nk)
    for ispin = 1:nspin
        for ik = 1:nk
            idx = (ispin - 1) * nk + ik
            band_offset = (ispin - 1) * minband
            for ib = 1:minband
                ebands[band_offset+ib, ik] = ebands_list[idx][ib]
            end
        end
    end

    # Check if mommat is available (all bands have velocity data)
    has_mommat = all(length(v) == 3 for vband in mommat_list for v in vband)

    mommat = nothing
    if has_mommat && !isempty(mommat_list[1][1])
        # Shape: (nspin * minband, 3, nk)
        mommat = zeros(nbands_total, 3, nk)
        for ispin = 1:nspin
            for ik = 1:nk
                idx = (ispin - 1) * nk + ik
                band_offset = (ispin - 1) * minband
                for ib = 1:minband
                    mommat[band_offset+ib, :, ik] = mommat_list[idx][ib]
                end
            end
        end
        # Convert Ry → Ha
        mommat .*= 0.5
    end

    # Convert Ry → Ha
    efermi = efermi_ry * 0.5
    ebands .*= 0.5

    # DOS weight
    dosweight = nspin == 1 ? 2.0 : 1.0

    return (
        fermi = efermi,
        dosweight = dosweight,
        kpoints = kpoints,
        ebands = ebands,
        mommat = mommat,
    )
end

#=
    _detect_gene_files(directory)

Detect GENE files in directory.

Returns `NamedTuple` with basename if GENE files found, `nothing` otherwise.
=#
function _detect_gene_files(directory::String)
    for f in readdir(directory)
        if endswith(f, ".structure") && isfile(joinpath(directory, f))
            basename = splitext(f)[1]
            energy_file = joinpath(directory, basename * ".energy")
            if isfile(energy_file)
                return (basename = basename,)
            end
        end
    end
    return nothing
end

"""
    load_gene(directory::String) -> [`DFTData`](@ref)

Load GENE/Generic format calculation results.

# Required files

  - `case.structure`: Crystal structure (lattice + atoms)
  - `case.energy`: Band energies and Fermi level

# Notes

  - For spin-polarized calculations, bands are concatenated (nspin × nbands)
  - All spin-polarized data is stored as [`DFTData`](@ref) with dosweight=1.0
  - Energy units: Rydberg → Hartree (×0.5)
  - Positions in .structure are Cartesian (Bohr), converted to fractional

# Returns

  - [`DFTData`](@ref) containing crystal structure and band structure data.
"""
function load_gene(directory::String)
    detected = _detect_gene_files(directory)
    if isnothing(detected)
        error("No GENE files found in $directory (need .structure and .energy)")
    end

    basename = detected.basename
    struct_file = joinpath(directory, basename * ".structure")
    energy_file = joinpath(directory, basename * ".energy")

    # Read structure
    struct_data = _read_gene_struct(struct_file)

    # Read energy
    energy_data = _read_gene_energy(energy_file)

    # Get number of electrons (not available in GENE format, estimate from nelect)
    # GENE doesn't store nelect, use -1 as placeholder
    nelect = -1.0

    # Reshape ebands to 3D: (nbands, nkpts, 1)
    nbands, nkpts = size(energy_data.ebands)
    ebands_3d = reshape(energy_data.ebands, nbands, nkpts, 1)

    # Weights: uniform (GENE doesn't store weights)
    weights = ones(nkpts) / nkpts

    # Occupations: zeros (GENE doesn't store occupations)
    occupations_3d = zeros(nbands, nkpts, 1)

    return DFTData(
        lattice = struct_data.lattice,
        positions = struct_data.positions,
        species = struct_data.species,
        kpoints = energy_data.kpoints,
        weights = weights,
        ebands = ebands_3d,
        occupations = occupations_3d,
        fermi = energy_data.fermi,
        nelect = nelect,
        magmom = nothing,
    )
end

#=
    _detect_gene(directory)

Format detection function for auto-loader.
Returns `true` if directory contains GENE files.
=#
function _detect_gene(directory::String)
    return !isnothing(_detect_gene_files(directory))
end
