# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

# Wien2k file loader for BoltzTraP.jl
# Reference: Python BoltzTraP2 dft.py (Wien2kLoader), io.py (W2Kene, W2Kfermi)
# Reference: ASE ase/io/wien2k.py (read_struct, c2p, coorsys)

using LinearAlgebra: I, inv

#=
    _ffloat(s)

Parse a Fortran-style floating point string.

Handles the "D" notation used in Fortran (e.g., "1.23D-05" → 1.23e-05).
=#
function _ffloat(s::AbstractString)
    parse(Float64, replace(lowercase(s), "d" => "e"))
end

#=
    _wien2k_c2p(lattice_type)

Return the conventional-to-primitive transformation matrix for Wien2k lattice types.

Based on ASE's ase.io.wien2k.c2p function.
Apply as: primitive_cell = c2p * conventional_cell
=#
function _wien2k_c2p(lattice_type::AbstractString)
    if lattice_type == "P"
        return Matrix{Float64}(I, 3, 3)
    elseif lattice_type == "F"
        return [0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]
    elseif lattice_type == "I"
        return [-0.5 0.5 0.5; 0.5 -0.5 0.5; 0.5 0.5 -0.5]
    elseif lattice_type == "C"
        return [0.5 0.5 0.0; 0.5 -0.5 0.0; 0.0 0.0 -1.0]
    elseif lattice_type == "B"
        return [0.5 0.0 0.5; 0.0 1.0 0.0; 0.5 0.0 -0.5]
    elseif lattice_type == "A"
        return [-1.0 0.0 0.0; 0.0 -0.5 0.5; 0.0 0.5 0.5]
    elseif lattice_type == "R"
        return [
            2.0/3.0 1.0/3.0 1.0/3.0;
            -1.0/3.0 1.0/3.0 1.0/3.0;
            -1.0/3.0 -2.0/3.0 1.0/3.0
        ]
    else
        error("Unknown Wien2k lattice type: $lattice_type")
    end
end

#=
    _wien2k_coorsys(a, b, c, alpha, beta, gamma)

Convert lattice constants to Cartesian cell vectors.

Based on ASE's ase.io.wien2k.coorsys function.
Returns 3x3 matrix with lattice vectors as columns.
=#
function _wien2k_coorsys(a::Real, b::Real, c::Real, alpha::Real, beta::Real, gamma::Real)
    # Convert angles from degrees to radians
    cal = cosd(alpha)
    cbe = cosd(beta)
    cga = cosd(gamma)
    sga = sind(gamma)

    # Construct cell matrix (ASE convention: rows are vectors, we transpose for columns)
    cell = [
        a b*cga c*cbe;
        0.0 b*sga c*(cal - cbe*cga)/sga;
        0.0 0.0 c*sqrt(1 - cal^2 - cbe^2 - cga^2 + 2*cal*cbe*cga)/sga
    ]
    return cell  # 3x3 with columns as lattice vectors
end

#=
    _read_wien2k_struct(filename)

Read Wien2k .struct file.

Returns named tuple with:
- lattice: 3x3 lattice vectors in Bohr (columns are vectors)
- positions: 3xN fractional coordinates
- species: Vector of element symbols
- lattice_type: Single character lattice type (P, F, I, etc.)
=#
function _read_wien2k_struct(filename::String)
    lines = readlines(filename)

    # Line 1: Title (ignored)
    # Line 2: Lattice type [0:3], Number of atoms [27:30] (0-indexed in Python)
    # In Julia 1-indexed: [1:3] and [28:30]
    lattice_raw = lines[2][1:3]  # Keep spaces for comparison
    nat = parse(Int, strip(lines[2][28:30]))

    # Convert Wien2k lattice type to standard notation
    # Compare with exact 3-char strings including spaces
    lattice_type = if startswith(lattice_raw, "P")
        "P"
    elseif startswith(lattice_raw, "H")
        "P"  # Hexagonal uses P with special angles
    elseif startswith(lattice_raw, "R")
        "R"
    elseif startswith(lattice_raw, "F")
        "F"
    elseif startswith(lattice_raw, "B")
        "I"  # Body-centered
    elseif lattice_raw == "CXY"
        "C"
    elseif lattice_raw == "CXZ"
        "B"
    elseif lattice_raw == "CYZ"
        "A"
    else
        error("Unknown Wien2k lattice type: '$lattice_raw'")
    end

    # Line 4 (index 4): Cell parameters (6 values, each 10 chars)
    # a, b, c in Bohr; alpha, beta, gamma in degrees
    cell_line = lines[4]
    a = parse(Float64, cell_line[1:10])
    b = parse(Float64, cell_line[11:20])
    c = parse(Float64, cell_line[21:30])
    alpha = parse(Float64, cell_line[31:40])
    beta = parse(Float64, cell_line[41:50])
    gamma = parse(Float64, cell_line[51:60])

    # Override angles for hexagonal lattice
    if startswith(lattice_raw, "H")
        alpha, beta, gamma = 90.0, 90.0, 120.0
    end

    # Construct conventional cell
    conv_cell = _wien2k_coorsys(a, b, c, alpha, beta, gamma)

    # Transform to primitive cell
    # For column-vector convention in Julia: prim_cell = conv_cell * c2p'
    # This is equivalent to Python's row-vector convention: prim_cell = c2p @ conv_cell
    c2p = _wien2k_c2p(lattice_type)
    prim_cell = conv_cell * transpose(c2p)

    # Parse atomic positions
    positions = Vector{Vector{Float64}}()
    species = String[]
    iline = 5  # Start after cell parameters (1-indexed)

    for _ = 1:nat
        # Position line: X [12:22], Y [25:35], Z [38:48] (0-indexed)
        # In Julia 1-indexed: X [13:22], Y [26:35], Z [39:48]
        pos_line = lines[iline]
        x = parse(Float64, pos_line[13:22])
        y = parse(Float64, pos_line[26:35])
        z = parse(Float64, pos_line[39:48])
        push!(positions, [x, y, z])
        iline += 1

        # MULT line: multiplicity [15:17] (0-indexed -> [16:17] in Julia)
        mult_line = lines[iline]
        mult = parse(Int, strip(mult_line[16:17]))
        iline += 1

        # Read equivalent positions (mult-1 additional positions)
        for _ = 2:mult
            pos_line = lines[iline]
            x = parse(Float64, pos_line[13:22])
            y = parse(Float64, pos_line[26:35])
            z = parse(Float64, pos_line[39:48])
            push!(positions, [x, y, z])
            iline += 1
        end

        # Species line: element symbol [1:2]
        species_line = lines[iline]
        element = strip(species_line[1:2])

        # Add species for all equivalent positions
        for _ = 1:mult
            push!(species, element)
        end

        # Skip remaining lines for this atom (3 more lines: NPT, LOCAL ROT MATRIX x3)
        iline += 4
    end

    # Convert positions to matrix (3 x natoms)
    pos_matrix = hcat(positions...)

    # Transform positions from conventional to primitive fractional coordinates
    # When the cell transforms as: prim_cell = c2p * conv_cell
    # Cartesian coordinates are preserved, so:
    #   cart = conv_cell * pos_conv = prim_cell * pos_prim
    #   pos_prim = inv(prim_cell) * conv_cell * pos_conv
    #            = inv(c2p * conv_cell) * conv_cell * pos_conv
    #            = inv(c2p) * pos_conv
    # But for rhombohedral (R), ASE scales atoms with the cell (different behavior)
    if lattice_type != "R"
        c2p_inv = inv(c2p)
        pos_matrix = c2p_inv * pos_matrix
    end

    # Wrap positions to [0, 1) for fractional coordinates (like ASE's atoms.wrap())
    pos_matrix = mod.(pos_matrix, 1.0)

    return (
        lattice = prim_cell,
        positions = pos_matrix,
        species = species,
        lattice_type = lattice_type,
    )
end

#=
    _read_wien2k_energy(filename, conv)

Read Wien2k .energy file.

Parameters:
- filename: Path to .energy or .energyso file
- conv: c2p transformation matrix from _wien2k_c2p

Returns named tuple with:
- kpoints: 3xN k-points in fractional coordinates (after c2p transformation)
- ebands: (nbands, nkpts) band energies in Hartree
=#
function _read_wien2k_energy(filename::String, conv::AbstractMatrix)
    lines = readlines(filename)

    # Find first valid k-point line
    linenumber = 1
    local nband_first::Int
    while linenumber <= length(lines)
        ll = lines[linenumber]
        try
            # K-point line: kx [1:19], ky [20:38], kz [39:57], nband [74:79]
            parse(Float64, ll[1:19])
            parse(Float64, ll[20:38])
            parse(Float64, ll[39:57])
            nband_first = parse(Int, strip(ll[74:79]))
            break
        catch
            linenumber += 1
        end
    end

    # Read all k-points and energies
    kpoints = Vector{Vector{Float64}}()
    ebands_list = Vector{Vector{Float64}}()
    minband = typemax(Int)

    while linenumber <= length(lines)
        ll = lines[linenumber]
        try
            # Parse k-point line
            kx = parse(Float64, ll[1:19])
            ky = parse(Float64, ll[20:38])
            kz = parse(Float64, ll[39:57])
            nband = parse(Int, strip(ll[74:79]))

            if nband < minband
                minband = nband
            end
            linenumber += 1

            # Read band energies
            eband = Float64[]
            for _ = 1:nband
                # Energy line: band_index energy (Rydberg)
                # Note: Wien2k uses Fortran "D" notation (e.g., "1.23D-05")
                parts = split(lines[linenumber])
                e = _ffloat(parts[2])
                push!(eband, e)
                linenumber += 1
            end

            push!(kpoints, [kx, ky, kz])
            push!(ebands_list, eband)
        catch
            break
        end
    end

    # Transform k-points: kpoints @ conv (Python) = conv' * kpoints (Julia)
    # Since conv is 3x3 and kpoints are column vectors, we use kpoints' * conv'
    kpts_matrix = hcat(kpoints...)  # 3 x nkpts
    kpts_transformed = conv' * kpts_matrix  # Apply transformation

    # Truncate all bands to minimum band count
    nkpts = length(ebands_list)
    ebands = zeros(minband, nkpts)
    for ik = 1:nkpts
        ebands[:, ik] = ebands_list[ik][1:minband]
    end

    # Convert from Rydberg to Hartree (multiply by 0.5)
    ebands .*= 0.5

    return (kpoints = kpts_transformed, ebands = ebands)
end

#=
    _read_wien2k_fermi(filename)

Read Fermi level from Wien2k .scf file.

Returns Fermi energy in Hartree.
=#
function _read_wien2k_fermi(filename::String)
    fermi = nothing
    for line in eachline(filename)
        if startswith(line, ":FER")
            # Fermi level is after the '=' sign
            # Format: ":FER  : F E R M I - ENERGY(...)=   value"
            # Python: 0.5 * float(l[38:53]) - but line length varies
            # Julia equivalent: start from index 39, read to end
            fermi_str = strip(line[39:end])
            fermi = 0.5 * parse(Float64, fermi_str)  # Rydberg to Hartree
        end
    end
    if isnothing(fermi)
        error("Fermi level not found in $filename")
    end
    return fermi
end

#=
    _read_wien2k_mommat_bounds(filename, nkpts)

Read band bounds (nemin, nemax) from Wien2k .mommat2 file.

This extracts the minimum and maximum band indices for which momentum matrix
elements are available. Used to trim ebands to match Python BoltzTraP2.

Returns (nemin, nemax) as 1-indexed integers (converted from Fortran indexing).
=#
function _read_wien2k_mommat_bounds(filename::String, nkpts::Int)
    lines = readlines(filename)

    # Skip first 2 header lines
    il = 3  # 1-indexed

    brk = Vector{Tuple{Int,Int}}()

    for _ = 1:nkpts
        # Parse nemin, nemax from header line
        # Format: "   KP:   N NEMIN NEMAX :  nemin nemax dE: ..."
        parts = split(lines[il])
        nemin_k = parse(Int, parts[6])
        nemax_k = parse(Int, parts[7])
        push!(brk, (nemin_k, nemax_k))

        # Skip to next k-point block
        # File format: header -> blank -> band_data (triangular) -> blank -> next header
        ne = nemax_k - nemin_k + 1
        # Number of band data lines: triangular format = ne*(ne+1)/2
        nlines_band = div(ne * (ne + 1), 2)
        # Advance: +1 (blank after header) + band_data + +1 (blank before next header)
        il += 1 + nlines_band + 1 + 1  # total = 2 + nlines_band + 1
    end

    # Find global nemin (max of local nemin) and nemax (min of local nemax)
    nemin = maximum(first.(brk))
    nemax = minimum(last.(brk))

    return (nemin, nemax)
end

#=
    _detect_wien2k_files(directory)

Detect Wien2k files in directory and return case name and energy file info.

Returns named tuple with:
- basename: Case name (e.g., "Si" for Si.struct)
- struct_file: Path to .struct file
- energy_file: Path to .energy or .energyso file
- scf_file: Path to .scf file
- mommat_file: Path to .mommat2 file or nothing
- dosweight: 2.0 for .energy, 1.0 for .energyso (SOC)
- nspin: 1 for non-spin-polarized, 2 for spin-polarized
=#
function _detect_wien2k_files(directory::String)
    # Find .struct file
    struct_files = filter(f -> endswith(f, ".struct"), readdir(directory))
    if isempty(struct_files)
        error("No .struct file found in $directory")
    end

    basename = splitext(struct_files[1])[1]
    struct_file = joinpath(directory, struct_files[1])

    # Check for energy files (priority: energyup/dn > energyso > energy)
    energyup = joinpath(directory, basename * ".energyup")
    energydn = joinpath(directory, basename * ".energydn")
    energysoup = joinpath(directory, basename * ".energysoup")
    energysodn = joinpath(directory, basename * ".energysodn")
    energyso = joinpath(directory, basename * ".energyso")
    energy = joinpath(directory, basename * ".energy")

    local energy_file::String
    local dosweight::Float64
    local nspin::Int

    if isfile(energyup) && isfile(energydn)
        # Spin-polarized without SOC
        energy_file = energyup  # Will need to read both
        dosweight = 1.0
        nspin = 2
        error("Spin-polarized Wien2k (.energyup/.energydn) not yet implemented")
    elseif isfile(energysoup) && isfile(energysodn)
        # Spin-polarized with SOC
        energy_file = energysoup  # Will need to read both
        dosweight = 1.0
        nspin = 2
        error(
            "Spin-polarized Wien2k with SOC (.energysoup/.energysodn) not yet implemented",
        )
    elseif isfile(energyso)
        # Non-spin-polarized with SOC
        energy_file = energyso
        dosweight = 1.0  # SOC doubles bands, so dosweight=1
        nspin = 1  # Actually stored as nspin=2 in DFTData due to dosweight
    elseif isfile(energy)
        # Non-spin-polarized without SOC
        energy_file = energy
        dosweight = 2.0
        nspin = 1
    else
        error("No .energy or .energyso file found for case $basename in $directory")
    end

    # Check for .scf file
    scf_file = joinpath(directory, basename * ".scf")
    if !isfile(scf_file)
        error("No .scf file found: $scf_file")
    end

    # Check for optional .mommat2 file
    mommat_file = joinpath(directory, basename * ".mommat2")
    mommat_file = isfile(mommat_file) ? mommat_file : nothing

    return (
        basename = basename,
        struct_file = struct_file,
        energy_file = energy_file,
        scf_file = scf_file,
        mommat_file = mommat_file,
        dosweight = dosweight,
        nspin = nspin,
    )
end

"""
    load_wien2k(directory) -> [`DFTData`](@ref)

Load Wien2k calculation data from directory.

# Required files

  - `case.struct`: Crystal structure
  - `case.energy` or `case.energyso`: Band energies
  - `case.scf`: Fermi level

# Notes

  - Spin-polarized calculations (.energyup/.energydn) are detected but not yet supported
  - SOC calculations (.energyso) use dosweight=1.0, resulting in NSpin=2 in [`DFTData`](@ref)
  - All energies converted from Rydberg to Hartree
  - Lattice vectors in Bohr

# Returns

  - [`DFTData`](@ref) where NSpin is determined by dosweight (2.0 → 1, 1.0 → 2)    # Detect files
"""
function load_wien2k(directory::String)
    # Detect files
    files = _detect_wien2k_files(directory)

    # Read structure
    struct_data = _read_wien2k_struct(files.struct_file)

    # Get c2p transformation matrix for k-point transformation
    # Note: BoltzTraP2 uses "P" (identity) for rhombohedral lattices
    lattice_type_for_kpts = struct_data.lattice_type == "R" ? "P" : struct_data.lattice_type
    c2p = _wien2k_c2p(lattice_type_for_kpts)

    # Read energies
    energy_data = _read_wien2k_energy(files.energy_file, c2p)

    # Read Fermi level
    fermi = _read_wien2k_fermi(files.scf_file)

    # Get band energy matrix
    ebands = energy_data.ebands
    nbands_raw, nkpts = size(ebands)

    # If mommat2 file exists, trim bands to the range with available derivatives
    # This matches Python BoltzTraP2's behavior
    if !isnothing(files.mommat_file)
        nemin, nemax = _read_wien2k_mommat_bounds(files.mommat_file, nkpts)
        # Fortran 1-indexed, Julia is also 1-indexed
        ebands = ebands[nemin:nemax, :]
    end

    nbands = size(ebands, 1)

    # K-point weights (uniform for Wien2k - will be normalized)
    weights = ones(nkpts) / nkpts

    # Determine number of electrons from bands and Fermi level
    # Wien2k doesn't provide nelect directly, so compute from occupations
    # nelect = dosweight * sum(occupied_bands * weights)
    occupancy = ebands .< fermi  # (nbands, nkpts) Bool matrix
    nelect = round(files.dosweight * sum(occupancy .* weights'))

    # Note on SOC (spin-orbit coupling):
    # - Python BoltzTraP2 uses dosweight=1.0 for SOC (vs 2.0 for non-polarized)
    # - But the ebands array shape is still (nbands, nkpts), NOT doubled
    # - The dosweight affects transport calculations, not band storage
    # For DFTData{NSpin}, we store SOC as NSpin=1 (single spin channel)
    # The transport calculations need to handle dosweight separately if needed
    ebands_3d = reshape(ebands, nbands, nkpts, 1)
    occupations_3d = zeros(nbands, nkpts, 1)

    return DFTData(
        lattice = struct_data.lattice,
        positions = struct_data.positions,
        species = struct_data.species,
        kpoints = energy_data.kpoints,
        weights = weights,
        ebands = ebands_3d,
        occupations = occupations_3d,
        fermi = fermi,
        nelect = nelect,
        magmom = nothing,
    )
end

#=
    _detect_wien2k(directory)

Check if directory contains Wien2k calculation files.
=#
function _detect_wien2k(directory::String)
    if !isdir(directory)
        return false
    end

    # Look for *.struct file
    for f in readdir(directory)
        if endswith(f, ".struct") && isfile(joinpath(directory, f))
            basename = splitext(f)[1]
            scf = joinpath(directory, basename * ".scf")
            energy = joinpath(directory, basename * ".energy")
            energyso = joinpath(directory, basename * ".energyso")
            if isfile(scf) && (isfile(energy) || isfile(energyso))
                return true
            end
        end
    end
    return false
end
