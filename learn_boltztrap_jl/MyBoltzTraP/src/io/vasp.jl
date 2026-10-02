# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2017-2025 Georg K. H. Madsen, Jesús Carrete, Matthieu J. Verstraete
# Copyright (C) 2026 Hiroharu Sugawara (Julia port)
# Part of BoltzTraP.jl - Port of BoltzTraP2/dft.py

# Unit conversion constants: EV_TO_HA, ANG_TO_BOHR defined in units.jl
# VASP uses eV and Å, internal uses Hartree and Bohr

#=
    read_poscar(filename)

Read VASP POSCAR/CONTCAR file.

# Returns
Named tuple with:
- `lattice`: 3×3 lattice vectors (columns)
- `species`: Element symbols
- `counts`: Number of atoms per species
- `positions`: 3×natoms atomic positions (columns)
- `comment`: Comment line
=#
function read_poscar(filename::String)
    lines = readlines(filename)

    # Line 1: Comment
    comment = lines[1]

    # Line 2: Scaling factor
    scale = parse(Float64, lines[2])

    # Lines 3-5: Lattice vectors (as columns)
    lattice = zeros(3, 3)
    for i = 1:3
        lattice[:, i] = parse.(Float64, split(lines[2+i]))
    end
    lattice .*= scale

    # Line 6: Element names (VASP 5+) or atom counts
    line6 = split(lines[6])
    if all(x -> !isnothing(tryparse(Int, x)), line6)
        # VASP 4 format: no element names
        species = nothing
        counts = parse.(Int, line6)
        pos_start = 7
    else
        # VASP 5+ format: element names present
        species = String.(line6)
        counts = parse.(Int, split(lines[7]))
        pos_start = 8
    end

    # Selective dynamics / coordinate mode
    mode_line = lines[pos_start]
    if startswith(lowercase(mode_line), "s")  # Selective dynamics
        pos_start += 1
        mode_line = lines[pos_start]
    end
    cartesian =
        startswith(lowercase(mode_line), "c") || startswith(lowercase(mode_line), "k")
    pos_start += 1

    # Atomic positions (as columns)
    natoms = sum(counts)
    positions = zeros(3, natoms)
    for i = 1:natoms
        coords = split(lines[pos_start+i-1])
        positions[:, i] = parse.(Float64, coords[1:3])
    end

    # Convert Direct (fractional) → Cartesian if needed
    if !cartesian
        positions = lattice * positions
    end

    return (
        lattice = lattice,
        species = species,
        counts = counts,
        positions = positions,
        comment = comment,
    )
end

#=
    read_vasprun(filename)

Read band structure data from VASP vasprun.xml file.

# Returns
Named tuple with:
- `lattice`: Lattice vectors (3, 3) in Ångström (columns are vectors)
- `rec_lattice`: Reciprocal lattice vectors (3, 3) in 1/Ångström
- `positions`: Atomic positions (3, natoms) in fractional coordinates
- `species`: Atom species names
- `kpoints`: K-points (3, nkpts) in fractional coordinates
- `weights`: K-point weights (nkpts,)
- `ebands`: Band energies (nbands, nkpts, nspin) in eV
- `occupations`: Occupations (nbands, nkpts, nspin)
- `fermi`: Fermi energy in eV
- `nelect`: Number of electrons
=#
function read_vasprun(filename::String)
    content = read(filename, String)

    # Extract Fermi energy (eV)
    fermi = 0.0
    m = match(r"<i name=\"efermi\">\s*([-\d.Ee+]+)\s*</i>", content)
    if !isnothing(m)
        fermi = parse(Float64, m.captures[1])
    end

    # Extract number of electrons
    nelect = 0.0
    m = match(r"<i name=\"NELECT\"[^>]*>\s*([-\d.Ee+]+)\s*</i>", content)
    if !isnothing(m)
        nelect = parse(Float64, m.captures[1])
    end

    # Extract lattice vectors from final structure
    # Look for the last <structure> block (finalpos)
    lattice = zeros(3, 3)
    structure_pattern = r"<structure[^>]*name=\"finalpos\"[^>]*>(.*?)</structure>"s
    m = match(structure_pattern, content)
    if isnothing(m)
        # Fall back to any structure block
        structure_pattern = r"<structure[^>]*>(.*?)</structure>"s
        m = match(structure_pattern, content)
    end

    if !isnothing(m)
        struct_content = m.captures[1]

        # Parse basis vectors
        basis_pattern = r"<varray name=\"basis\"[^>]*>(.*?)</varray>"s
        bm = match(basis_pattern, struct_content)
        if !isnothing(bm)
            row_pattern = r"<v>(.*?)</v>"
            for (i, rm) in enumerate(eachmatch(row_pattern, bm.captures[1]))
                if i <= 3
                    vals = parse.(Float64, split(strip(rm.captures[1])))
                    lattice[:, i] = vals  # Store as columns
                end
            end
        end
    end

    # Extract reciprocal lattice vectors
    rec_lattice = zeros(3, 3)
    rec_pattern = r"<varray name=\"rec_basis\"[^>]*>(.*?)</varray>"s
    m = match(rec_pattern, content)
    if !isnothing(m)
        row_pattern = r"<v>(.*?)</v>"
        for (i, rm) in enumerate(eachmatch(row_pattern, m.captures[1]))
            if i <= 3
                vals = parse.(Float64, split(strip(rm.captures[1])))
                rec_lattice[:, i] = vals
            end
        end
    end

    # Extract atomic positions (fractional)
    positions = zeros(3, 0)
    pos_pattern = r"<varray name=\"positions\"[^>]*>(.*?)</varray>"s
    # Find the last occurrence (finalpos)
    all_pos = collect(eachmatch(pos_pattern, content))
    if !isempty(all_pos)
        pos_content = all_pos[end].captures[1]
        row_pattern = r"<v>(.*?)</v>"
        for rm in eachmatch(row_pattern, pos_content)
            vals = parse.(Float64, split(strip(rm.captures[1])))
            positions = hcat(positions, vals)
        end
    end

    # Extract atom species
    species = String[]
    atoms_pattern = r"<array name=\"atoms\"[^>]*>.*?<set>(.*?)</set>"s
    m = match(atoms_pattern, content)
    if !isnothing(m)
        # Handle both compact (<rc><c>X</c>...) and pretty-printed (<rc>\n<c>X</c>...) formats
        rc_pattern = r"<rc>\s*<c>\s*(\S+)\s*</c>"s
        for rm in eachmatch(rc_pattern, m.captures[1])
            push!(species, strip(rm.captures[1]))
        end
    end

    # Extract k-points
    kpoints = zeros(3, 0)
    kpt_pattern = r"<varray name=\"kpointlist\"[^>]*>(.*?)</varray>"s
    m = match(kpt_pattern, content)
    if !isnothing(m)
        row_pattern = r"<v>(.*?)</v>"
        for rm in eachmatch(row_pattern, m.captures[1])
            vals = parse.(Float64, split(strip(rm.captures[1])))
            kpoints = hcat(kpoints, vals)
        end
    end
    nkpts = size(kpoints, 2)

    # Extract k-point weights
    weights = zeros(nkpts)
    weights_pattern = r"<varray name=\"weights\"[^>]*>(.*?)</varray>"s
    m = match(weights_pattern, content)
    if !isnothing(m)
        row_pattern = r"<v>(.*?)</v>"
        for (i, rm) in enumerate(eachmatch(row_pattern, m.captures[1]))
            if i <= nkpts
                weights[i] = parse(Float64, strip(rm.captures[1]))
            end
        end
    end

    # Extract eigenvalues and occupations
    # Structure: <eigenvalues><array><set><set spin="1"><set>...</set></set></set></array></eigenvalues>
    eigen_pattern =
        r"<eigenvalues>\s*<array>.*?<set>(.*?)</set>\s*</array>\s*</eigenvalues>"s
    m = match(eigen_pattern, content)

    ebands = Array{Float64,3}(undef, 0, 0, 0)
    occupations = Array{Float64,3}(undef, 0, 0, 0)

    if !isnothing(m)
        eigenset_content = m.captures[1]

        # Find spin sets
        spin_pattern =
            r"<set comment=\"spin\s*(\d+)\"[^>]*>(.*?)</set>\s*(?=<set comment=\"spin|$)"s
        spin_matches = collect(eachmatch(spin_pattern, eigenset_content))

        if isempty(spin_matches)
            # No spin polarization - single spin channel
            spin_matches = [(captures = ["1", eigenset_content],)]
        end

        nspin = length(spin_matches)

        # Parse first k-point to get nbands
        first_kpt_pattern = r"<set comment=\"kpoint\s*\d+\"[^>]*>(.*?)</set>"s
        first_kpt = match(first_kpt_pattern, spin_matches[1].captures[2])
        if !isnothing(first_kpt)
            r_pattern = r"<r>(.*?)</r>"
            nbands = length(collect(eachmatch(r_pattern, first_kpt.captures[1])))

            ebands = zeros(nbands, nkpts, nspin)
            occupations = zeros(nbands, nkpts, nspin)

            for (ispin, spin_match) in enumerate(spin_matches)
                spin_content = spin_match.captures[2]

                # Parse each k-point
                kpt_pattern = r"<set comment=\"kpoint\s*(\d+)\"[^>]*>(.*?)</set>"s
                for km in eachmatch(kpt_pattern, spin_content)
                    ik = parse(Int, km.captures[1])
                    kpt_content = km.captures[2]

                    # Parse energy and occupation pairs
                    r_pattern = r"<r>\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*</r>"
                    for (ib, rm) in enumerate(eachmatch(r_pattern, kpt_content))
                        if ib <= nbands && ik <= nkpts
                            ebands[ib, ik, ispin] = parse(Float64, rm.captures[1])
                            occupations[ib, ik, ispin] = parse(Float64, rm.captures[2])
                        end
                    end
                end
            end
        end
    end

    # Extract MAGMOM (magnetic moments)
    # Format in vasprun.xml: <v name="MAGMOM">  1.0  1.0  ... </v>
    magmom = nothing
    magmom_pattern = r"<v name=\"MAGMOM\">\s*(.*?)\s*</v>"s
    m = match(magmom_pattern, content)
    if !isnothing(m)
        magmom_str = strip(m.captures[1])
        if !isempty(magmom_str)
            magmom = parse.(Float64, split(magmom_str))
        end
    end

    return (
        lattice = lattice,
        rec_lattice = rec_lattice,
        positions = positions,
        species = species,
        kpoints = kpoints,
        weights = weights,
        ebands = ebands,
        occupations = occupations,
        fermi = fermi,
        nelect = nelect,
        magmom = magmom,
    )
end

"""
    load_vasp(directory) -> [`DFTData`](@ref)

Load VASP calculation data from directory and convert to atomic units.

Reads vasprun.xml for band structure and structure information.
Converts from VASP units (eV, Å) to atomic units (Hartree, Bohr).

# Returns

  - [`DFTData`](@ref) where NSpin is 1 (non-spin-polarized) or 2 (spin-polarized):
  - `lattice`: Lattice vectors (3, 3) in Bohr (columns are vectors)
  - `positions`: Atomic positions (3, natoms) in fractional coordinates
  - `species`: Atom species names
  - `kpoints`: K-points (3, nkpts) in fractional coordinates
  - `weights`: K-point weights (nkpts,)
  - `ebands`: Band energies (nbands, nkpts, nspin) in Hartree
  - `occupations`: Occupations (nbands, nkpts, nspin)
  - `fermi`: Fermi energy in Hartree
  - `nelect`: Number of electrons
  - `magmom`: Magnetic moments (nothing for non-spin-polarized)
"""
function load_vasp(directory::String)
    vasprun = joinpath(directory, "vasprun.xml")
    if !isfile(vasprun)
        error("vasprun.xml not found in $directory")
    end
    raw = read_vasprun(vasprun)

    # Convert to atomic units (Hartree, Bohr)
    return DFTData(
        lattice = raw.lattice * ANG_TO_BOHR,
        positions = raw.positions,  # fractional, no conversion needed
        species = raw.species,
        kpoints = raw.kpoints,      # fractional, no conversion needed
        weights = raw.weights,
        ebands = raw.ebands * EV_TO_HA,
        occupations = raw.occupations,
        fermi = raw.fermi * EV_TO_HA,
        nelect = raw.nelect,
        magmom = raw.magmom,  # Extracted from vasprun.xml
    )
end

#=
    read_eigenval(filename)

Read VASP EIGENVAL file.

# Returns
Named tuple with:
- `nelect`: Number of electrons
- `nkpts`: Number of k-points
- `nbands`: Number of bands
- `kpoints`: K-points (3, nkpts) in fractional coordinates
- `weights`: K-point weights (nkpts,)
- `ebands`: Band energies (nbands, nkpts, nspin) in eV
- `occupations`: Occupations (nbands, nkpts, nspin)
=#
function read_eigenval(filename::String)
    lines = readlines(filename)

    # Line 1: system info (natoms, ?, ?, ispin)
    header1 = parse.(Int, split(strip(lines[1])))
    nspin = header1[end]  # 1 or 2

    # Lines 2-5: additional info (skip)
    # Line 6: NELECT, NKPTS, NBANDS
    header6 = split(strip(lines[6]))
    nelect = parse(Int, header6[1])
    nkpts = parse(Int, header6[2])
    nbands = parse(Int, header6[3])

    kpoints = zeros(3, nkpts)
    weights = zeros(nkpts)
    ebands = zeros(nbands, nkpts, nspin)
    occupations = zeros(nbands, nkpts, nspin)

    line_idx = 8  # Skip header and blank line

    for ik = 1:nkpts
        # K-point line: kx ky kz weight
        kline = split(strip(lines[line_idx]))
        kpoints[1, ik] = parse(Float64, kline[1])
        kpoints[2, ik] = parse(Float64, kline[2])
        kpoints[3, ik] = parse(Float64, kline[3])
        weights[ik] = parse(Float64, kline[4])
        line_idx += 1

        # Band lines: band_index energy [occupation] [energy_spin2] [occupation_spin2]
        for ib = 1:nbands
            bline = split(strip(lines[line_idx]))
            if nspin == 1
                # Format: band energy occupation
                ebands[ib, ik, 1] = parse(Float64, bline[2])
                if length(bline) >= 3
                    occupations[ib, ik, 1] = parse(Float64, bline[3])
                end
            else
                # Format: band energy_up energy_down
                ebands[ib, ik, 1] = parse(Float64, bline[2])
                ebands[ib, ik, 2] = parse(Float64, bline[3])
            end
            line_idx += 1
        end

        line_idx += 1  # Skip blank line between k-points
    end

    return (
        nelect = nelect,
        nkpts = nkpts,
        nbands = nbands,
        kpoints = kpoints,
        weights = weights,
        ebands = ebands,
        occupations = occupations,
    )
end
