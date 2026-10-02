# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

using LinearAlgebra: I, det

# Chemical symbols for element name validation
const CHEMICAL_SYMBOLS = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

#=
    _normalize_element_case(name)

Normalize element name to proper case (first letter uppercase, rest lowercase).
=#
function _normalize_element_case(name::AbstractString)
    if isempty(name)
        return ""
    end
    return uppercase(name[1:1]) * lowercase(name[2:end])
end

#=
    _unpack_qe_element(name)

Extract chemical element symbol from QE extended element name.

QE allows element names like "Fe1", "Cr_up", "C-h" etc.
Returns the base element symbol (e.g., "Fe", "Cr", "C").
=#
function _unpack_qe_element(name::AbstractString)
    name = strip(name)
    if length(name) > 3
        error("Quantum ESPRESSO element name $name is too long")
    end

    # Check for dash or underscore separator
    if occursin('_', name) || occursin('-', name)
        fields = split(replace(name, '-' => '_'), '_')
        if length(fields) != 2
            error("Invalid Quantum ESPRESSO element name $name")
        end
        element_name = _normalize_element_case(fields[1])
        if !(element_name in CHEMICAL_SYMBOLS)
            error("Invalid chemical element $element_name in QE name $name")
        end
        return element_name
    end

    # No separator: find largest valid element match from start
    normalized = _normalize_element_case(name)
    for i = length(normalized):-1:1
        substring = normalized[1:i]
        if substring in CHEMICAL_SYMBOLS
            return substring
        end
    end

    error("No valid chemical symbol found in QE element name $name")
end

#=
    _generate_qe_tags(extended_species)

Generate atomic tags from QE extended species names.

Python BoltzTraP2 uses these tags to generate dummy magmom values for symmetry breaking.
Same-named extended species get the same tag (0-based).

# Example
```julia
_generate_qe_tags(["Fe1", "Fe1", "Fe2", "Fe2"]) # → [0, 0, 1, 1]
```

# Reference
Python BoltzTraP2 io.py:1413-1422 and io.py:1526-1528
=#
function _generate_qe_tags(extended_species::Vector{String})
    unique_ext = unique(extended_species)
    return [findfirst(==(s), unique_ext) - 1 for s in extended_species]  # 0-based
end

#=
    read_qe_bands_gnu(filename)

Read band energies from QE bands.x output (bands.dat.gnu format).

The GNU format has columns: k_distance, energy (eV)
Bands are separated by blank lines.

# Returns
- `kdist`: k-point distances along path (nkpts,)
- `ebands`: Band energies in eV (nbands, nkpts)
=#
function read_qe_bands_gnu(filename::String)
    lines = readlines(filename)

    # Parse data, grouping by bands (blank lines separate bands)
    bands_data = Vector{Vector{Tuple{Float64,Float64}}}()
    current_band = Vector{Tuple{Float64,Float64}}()

    for line in lines
        stripped = strip(line)
        if isempty(stripped)
            if !isempty(current_band)
                push!(bands_data, current_band)
                current_band = Vector{Tuple{Float64,Float64}}()
            end
        else
            parts = split(stripped)
            if length(parts) >= 2
                kdist = parse(Float64, parts[1])
                energy = parse(Float64, parts[2])
                push!(current_band, (kdist, energy))
            end
        end
    end
    # Don't forget the last band
    if !isempty(current_band)
        push!(bands_data, current_band)
    end

    # Convert to arrays
    nbands = length(bands_data)
    nkpts = length(bands_data[1])

    kdist = [bands_data[1][ik][1] for ik = 1:nkpts]
    ebands = zeros(nbands, nkpts)
    for ib = 1:nbands
        for ik = 1:nkpts
            ebands[ib, ik] = bands_data[ib][ik][2]
        end
    end

    return (kdist = kdist, ebands = ebands)
end

#=
    read_qe_kpoints(filename)

Read k-points from QE bands.x output (bands.dat format).

The bands.dat format header contains:
- Line 1: nbnd, nks
- Line 2+: k-point coordinates (crystal) followed by energies

# Returns
- `kpoints`: K-points in crystal coordinates (3, nkpts)
- `ebands`: Band energies in eV (nbands, nkpts)
=#
function read_qe_kpoints(filename::String)
    lines = readlines(filename)

    # Parse header
    header = split(strip(lines[1]))
    nbnd = parse(Int, header[1])
    nks = parse(Int, header[2])

    kpoints = zeros(3, nks)
    ebands = zeros(nbnd, nks)

    line_idx = 2
    for ik = 1:nks
        # K-point line: kx ky kz
        kline = split(strip(lines[line_idx]))
        kpoints[1, ik] = parse(Float64, kline[1])
        kpoints[2, ik] = parse(Float64, kline[2])
        kpoints[3, ik] = parse(Float64, kline[3])
        line_idx += 1

        # Energy lines (10 values per line typically)
        energies = Float64[]
        while length(energies) < nbnd
            eline = split(strip(lines[line_idx]))
            for e in eline
                push!(energies, parse(Float64, e))
            end
            line_idx += 1
        end
        ebands[:, ik] = energies[1:nbnd]
    end

    return (kpoints = kpoints, ebands = ebands)
end

#=
    read_qe_xml(filename)

Read structure and band data from QE data-file-schema.xml (QE 6.x+).

# Returns
Named tuple with:
- `lattice`: Lattice vectors (3, 3) in Bohr
- `positions`: Atomic positions (3, natoms) in fractional coordinates
- `species`: Atom species names
- `kpoints`: K-points (3, nkpts) in fractional coordinates
- `weights`: K-point weights (nkpts,)
- `ebands`: Band energies (nbands, nkpts, nspin) in Hartree
- `occupations`: Occupations (nbands, nkpts, nspin)
- `fermi`: Fermi energy in Hartree
- `nelect`: Number of electrons
- `nspin`: Number of spin channels (1 or 2)
=#
function read_qe_xml(filename::String)
    # Requires EzXML.jl for proper XML parsing
    # For now, provide a simple regex-based parser for essential data

    content = read(filename, String)

    # Extract lattice vectors (in Bohr) from output section
    lattice = zeros(3, 3)
    # Look for <cell> block in output/atomic_structure
    cell_match = match(r"<output>.*?<atomic_structure[^>]*>.*?<cell>(.*?)</cell>"s, content)
    if !isnothing(cell_match)
        cell_content = cell_match.captures[1]
        for i = 1:3
            pattern = Regex("<a$i>(.*?)</a$i>")
            m = match(pattern, cell_content)
            if !isnothing(m)
                vals = parse.(Float64, split(strip(m.captures[1])))
                lattice[:, i] = vals
            end
        end
    else
        # Fallback: use first occurrence
        for i = 1:3
            pattern = Regex("<a$i>(.*?)</a$i>")
            m = match(pattern, content)
            if !isnothing(m)
                vals = parse.(Float64, split(strip(m.captures[1])))
                lattice[:, i] = vals
            end
        end
    end

    # Extract reciprocal lattice vectors (in 2π/Bohr) for k-point transformation
    rlattvec = zeros(3, 3)
    rlat_match = match(r"<reciprocal_lattice>(.*?)</reciprocal_lattice>"s, content)
    if !isnothing(rlat_match)
        rlat_content = rlat_match.captures[1]
        for i = 1:3
            pattern = Regex("<b$i>(.*?)</b$i>")
            m = match(pattern, rlat_content)
            if !isnothing(m)
                vals = parse.(Float64, split(strip(m.captures[1])))
                rlattvec[:, i] = vals
            end
        end
    else
        # Compute reciprocal lattice from direct lattice if not present
        # b_i = 2π * (a_j × a_k) / (a_i · (a_j × a_k)) for cyclic i,j,k
        # Using the formula: reciprocal = 2π * inv(direct)'
        # But since we just need it for k-point transformation (which divides by 2π anyway),
        # we can use a simplified version
        if det(lattice) != 0
            # Standard formula: B = 2π * inv(A)'  but QE stores B differently
            # Actually QE's reciprocal lattice from XML is just inv(A)' without 2π factor
            # For the transformation, we use: k_frac = rlattvec \ k_cart
            # If rlattvec not provided, assume k-points are already fractional
            rlattvec = Matrix{Float64}(I, 3, 3)  # Identity means no transformation
        end
    end

    # Detect spin polarization from band_structure section
    spin_polarized = false
    lsda_match = match(r"<band_structure>.*?<lsda>(.*?)</lsda>"s, content)
    if !isnothing(lsda_match)
        spin_polarized = strip(lsda_match.captures[1]) == "true"
    end
    nspin = spin_polarized ? 2 : 1

    # Extract Fermi energy (in Hartree)
    # For semiconductors, QE may not output fermi_energy, use HOMO/LUMO average
    fermi = 0.0
    m = match(r"<fermi_energy>(.*?)</fermi_energy>", content)
    if !isnothing(m)
        fermi = parse(Float64, strip(m.captures[1]))
    else
        # Try HOMO/LUMO for semiconductors
        m_hol = match(r"<highestOccupiedLevel>(.*?)</highestOccupiedLevel>", content)
        m_lul = match(r"<lowestUnoccupiedLevel>(.*?)</lowestUnoccupiedLevel>", content)
        if !isnothing(m_hol) && !isnothing(m_lul)
            hol = parse(Float64, strip(m_hol.captures[1]))
            lul = parse(Float64, strip(m_lul.captures[1]))
            fermi = 0.5 * (hol + lul)
        elseif !isnothing(m_hol)
            fermi = parse(Float64, strip(m_hol.captures[1]))
        end
    end

    # Extract number of electrons
    nelect = 0.0
    m = match(r"<nelec>(.*?)</nelec>", content)
    if !isnothing(m)
        nelect = parse(Float64, strip(m.captures[1]))
    end

    # Extract number of k-points
    nk = 0
    m = match(r"<nks>(.*?)</nks>", content)
    if !isnothing(m)
        nk = parse(Int, strip(m.captures[1]))
    end

    # Extract number of bands (handle spin-polarized case)
    nbnd = 0
    if spin_polarized
        m = match(r"<nbnd_up>(.*?)</nbnd_up>", content)
        if !isnothing(m)
            nbnd = parse(Int, strip(m.captures[1]))
        end
    else
        m = match(r"<nbnd>(.*?)</nbnd>", content)
        if !isnothing(m)
            nbnd = parse(Int, strip(m.captures[1]))
        end
    end

    # Extract k-points (Cartesian) and eigenvalues
    kpoints_cart = zeros(3, nk)
    weights = zeros(nk)
    # For spin-polarized: Python BoltzTraP2 flattens spin into band dimension
    ebands_raw = zeros(nbnd * nspin, nk)
    occupations_raw = zeros(nbnd * nspin, nk)

    # Parse each ks_energies block
    ks_pattern = r"<ks_energies>(.*?)</ks_energies>"s
    for (ik, m) in enumerate(eachmatch(ks_pattern, content))
        block = m.captures[1]

        # K-point coordinates (Cartesian in units of 2π/a)
        km = match(r"<k_point[^>]*>(.*?)</k_point>", block)
        if !isnothing(km)
            kvals = parse.(Float64, split(strip(km.captures[1])))
            kpoints_cart[:, ik] = kvals
        end

        # Weight
        wm = match(r"weight=\"([^\"]+)\"", block)
        if !isnothing(wm)
            weights[ik] = parse(Float64, wm.captures[1])
        end

        # Eigenvalues
        em = match(r"<eigenvalues[^>]*>(.*?)</eigenvalues>"s, block)
        if !isnothing(em)
            evals = parse.(Float64, split(strip(em.captures[1])))
            ebands_raw[:, ik] = evals[1:(nbnd*nspin)]
        end

        # Occupations
        om = match(r"<occupations[^>]*>(.*?)</occupations>"s, block)
        if !isnothing(om)
            ovals = parse.(Float64, split(strip(om.captures[1])))
            occupations_raw[:, ik] = ovals[1:(nbnd*nspin)]
        end
    end

    # Transform k-points from Cartesian to fractional coordinates
    # k_frac = rlattvec \ k_cart (solve rlattvec * k_frac = k_cart)
    kpoints = zeros(3, nk)
    for ik = 1:nk
        kpoints[:, ik] = rlattvec \ kpoints_cart[:, ik]
    end
    # Wrap to [-0.5, 0.5] like Python BoltzTraP2
    kpoints .-= round.(kpoints)

    # Extract atomic positions from first <atomic_structure> block in output
    # (XML may contain multiple: input and output structures)
    positions_cart = zeros(3, 0)
    species = String[]
    extended_species = String[]  # Keep original extended names for tags generation

    # Look for atomic_positions in output section
    output_struct_match = match(
        r"<output>.*?<atomic_structure[^>]*>.*?<atomic_positions>(.*?)</atomic_positions>"s,
        content,
    )
    if !isnothing(output_struct_match)
        atom_pattern = r"<atom[^>]*name=\"([^\"]+)\"[^>]*>(.*?)</atom>"s
        for m in eachmatch(atom_pattern, output_struct_match.captures[1])
            raw_name = strip(m.captures[1])
            # Keep extended name for tags generation (Python BoltzTraP2 compatibility)
            push!(extended_species, raw_name)
            # Extract base element from QE extended name (e.g., "Cr1" -> "Cr")
            push!(species, _unpack_qe_element(raw_name))
            pos = parse.(Float64, split(strip(m.captures[2])))
            positions_cart = hcat(positions_cart, pos)
        end
    else
        # Fallback: first atomic_structure block
        atomic_struct_match = match(
            r"<atomic_structure[^>]*>.*?<atomic_positions>(.*?)</atomic_positions>"s,
            content,
        )
        if !isnothing(atomic_struct_match)
            atom_pattern = r"<atom[^>]*name=\"([^\"]+)\"[^>]*>(.*?)</atom>"s
            for m in eachmatch(atom_pattern, atomic_struct_match.captures[1])
                raw_name = strip(m.captures[1])
                # Keep extended name for tags generation (Python BoltzTraP2 compatibility)
                push!(extended_species, raw_name)
                # Extract base element from QE extended name (e.g., "Cr1" -> "Cr")
                push!(species, _unpack_qe_element(raw_name))
                pos = parse.(Float64, split(strip(m.captures[2])))
                positions_cart = hcat(positions_cart, pos)
            end
        end
    end

    # Convert atomic positions from Cartesian (Bohr) to fractional coordinates
    # ASE convention: pos_frac = pos_cart @ inv(cell) where cell has rows as vectors
    # In Julia column notation with lattice (columns as vectors):
    #   cell = lattice', so inv(cell') = inv(lattice)
    #   pos_frac = inv(lattice) @ pos_cart = lattice \ pos_cart
    # Then wrap to [0, 1) range
    natoms = size(positions_cart, 2)
    positions = zeros(3, natoms)
    for ia = 1:natoms
        frac = lattice \ positions_cart[:, ia]
        positions[:, ia] = frac .- floor.(frac)  # Wrap to [0, 1)
    end

    return (
        lattice = lattice,
        positions = positions,
        species = species,
        extended_species = extended_species,  # For tags generation (Python BoltzTraP2 compatibility)
        kpoints = kpoints,
        weights = weights,
        ebands = ebands_raw,
        occupations = occupations_raw,
        fermi = fermi,
        nelect = nelect,
        nspin = nspin,
    )
end

#=
    read_qe_output(filename)

Read band energies from QE pw.x standard output (stdout).

Parses the eigenvalue blocks printed by pw.x during nscf or bands calculations:
```
          k = 0.0000 0.0000 0.0000 (   283 PWs)   bands (ev):

    -5.8094   6.2539   6.2539   6.2539
```

# Returns
- `kpoints`: K-points in crystal coordinates (3, nkpts)
- `ebands`: Band energies in eV (nbands, nkpts)
=#
function read_qe_output(filename::String)
    content = read(filename, String)
    lines = split(content, '\n')

    kpoints_list = Vector{Vector{Float64}}()
    ebands_list = Vector{Vector{Float64}}()

    i = 1
    while i <= length(lines)
        line = lines[i]

        # Match k-point line: "k = 0.0000 0.0000 0.0000 (   283 PWs)   bands (ev):"
        m = match(r"k\s*=\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line)
        if !isnothing(m)
            # Extract k-point
            kpt = [parse(Float64, m.captures[j]) for j = 1:3]
            push!(kpoints_list, kpt)

            # Skip to next non-empty line (eigenvalues)
            i += 1
            energies = Float64[]

            # Read energy lines until next k-point or end of block
            while i <= length(lines)
                eline = strip(lines[i])

                # Stop if we hit another k-point or specific markers
                if startswith(eline, "k =") ||
                   startswith(eline, "the Fermi") ||
                   startswith(eline, "highest") ||
                   startswith(eline, "Writing") ||
                   startswith(eline, "End of") ||
                   occursin("CPU", eline)
                    break
                end

                # Parse energy values if line has numbers
                if !isempty(eline)
                    parts = split(eline)
                    for p in parts
                        try
                            push!(energies, parse(Float64, p))
                        catch
                            # Not a number, skip
                        end
                    end
                end
                i += 1
            end

            if !isempty(energies)
                push!(ebands_list, energies)
            end
        else
            i += 1
        end
    end

    # Convert to arrays
    nkpts = length(kpoints_list)
    if nkpts == 0
        error("No k-points found in QE output")
    end

    nbands = length(ebands_list[1])
    kpoints = zeros(3, nkpts)
    ebands = zeros(nbands, nkpts)

    for ik = 1:nkpts
        kpoints[:, ik] = kpoints_list[ik]
        ebands[:, ik] = ebands_list[ik][1:nbands]
    end

    return (kpoints = kpoints, ebands = ebands)
end

"""
    load_qe(directory) -> [`DFTData`](@ref)

Load Quantum ESPRESSO calculation data from directory.

Searches for XML data file (data-file-schema.xml or *.xml in out/ subdirectory).
Returns [`DFTData`](@ref) for non-magnetic or spin-polarized.
All data in atomic units (Hartree, Bohr).

# Arguments

  - `directory`: Path to QE output directory containing XML file

# Returns

  - [`DFTData`](@ref) with NSpin ∈ {1, 2}
"""
function load_qe(directory::String)
    # Find XML file
    xml_file = nothing

    # Try common locations
    candidates =
        [joinpath(directory, "data-file-schema.xml"), joinpath(directory, "data-file.xml")]

    # Search for *.xml files directly in the directory
    if isdir(directory)
        for f in readdir(directory)
            if endswith(f, ".xml")
                push!(candidates, joinpath(directory, f))
            end
        end
    end

    # Also search in out/ subdirectory
    out_dir = joinpath(directory, "out")
    if isdir(out_dir)
        for f in readdir(out_dir)
            if endswith(f, ".xml")
                push!(candidates, joinpath(out_dir, f))
            end
        end
    end

    # Find first existing file
    for candidate in candidates
        if isfile(candidate)
            xml_file = candidate
            break
        end
    end

    if isnothing(xml_file)
        error("No QE XML file found in $directory")
    end

    raw = read_qe_xml(xml_file)

    # QE XML already returns data in atomic units (Ha, Bohr)
    # ebands/occupations are already (nbands*nspin, nkpts) from read_qe_xml
    # Reshape to (nbands, nkpts, nspin) for compatibility with VASP format
    nspin_val = raw.nspin
    nbands_total = size(raw.ebands, 1)
    nbands = nbands_total ÷ nspin_val
    nkpts = size(raw.ebands, 2)

    # For spin-polarized: bands are stored as [up_bands..., down_bands...] per k-point
    # Need to reshape to (nbands, nkpts, nspin)
    ebands_3d = zeros(nbands, nkpts, nspin_val)
    occupations_3d = zeros(nbands, nkpts, nspin_val)
    for ispin = 1:nspin_val
        band_start = (ispin - 1) * nbands + 1
        band_end = ispin * nbands
        ebands_3d[:, :, ispin] = raw.ebands[band_start:band_end, :]
        occupations_3d[:, :, ispin] = raw.occupations[band_start:band_end, :]
    end

    # Generate dummy magmom from atomic tags (Python BoltzTraP2 compatibility)
    # Purpose: break symmetry for magnetic calculations
    # Reference: Python BoltzTraP2 io.py:1526-1528
    #   nruter["magmom"] = (np.array(nruter["atoms"].get_tags(), dtype=np.float64) + 1.0)
    magmom = nothing
    if !isempty(raw.extended_species)
        tags = _generate_qe_tags(raw.extended_species)
        magmom = Float64.(tags) .+ 1.0
    end

    # Return DFTData with appropriate type parameter
    return DFTData(
        lattice = raw.lattice,
        positions = raw.positions,
        species = raw.species,
        kpoints = raw.kpoints,
        weights = raw.weights,
        ebands = ebands_3d,
        occupations = occupations_3d,
        fermi = raw.fermi,
        nelect = Float64(raw.nelect),
        magmom = magmom,
    )
end
