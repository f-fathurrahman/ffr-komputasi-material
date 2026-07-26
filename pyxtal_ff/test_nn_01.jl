using AtomsBase, AtomsIO, NeighbourLists, Unitful, LinearAlgebra

# Atomic number dictionary (extend as needed)
const ATOMIC_NUMBER = Dict(
    "H" => 1,  "He"=>2,  "Li"=>3,  "Be"=>4,  "B"=>5,  "C"=>6,  "N"=>7,  "O"=>8,
    "F" => 9,  "Ne"=>10, "Na"=>11, "Mg"=>12, "Al"=>13, "Si"=>14, "P"=>15, "S"=>16,
    "Cl"=>17, "Ar"=>18, "K"=>19,  "Ca"=>20, "Sc"=>21, "Ti"=>22, "V"=>23,  "Cr"=>24,
    "Mn"=>25, "Fe"=>26, "Co"=>27, "Ni"=>28, "Cu"=>29, "Zn"=>30, "Ga"=>31, "Ge"=>32,
    "As"=>33, "Se"=>34, "Br"=>35, "Kr"=>36, "Rb"=>37, "Sr"=>38, "Y"=>39,  "Zr"=>40,
    "Nb"=>41, "Mo"=>42, "Tc"=>43, "Ru"=>44, "Rh"=>45, "Pd"=>46, "Ag"=>47, "Cd"=>48,
    "In"=>49, "Sn"=>50, "Sb"=>51, "Te"=>52, "I"=>53,  "Xe"=>54, "Cs"=>55, "Ba"=>56
)

"""
    read_first_frame(filename::String) -> AbstractSystem

Read the first structure from an extended XYZ file using AtomsIO.
"""
function read_first_frame(filename::String)
    # load_system returns an AbstractSystem with the first frame
    system = load_system(filename)
    return system
end

"""
    build_neighbor_data(system, rcut; weight_on=false)

Build all neighbour-related arrays used in the original Python script.

Returns a tuple:
    (center_atoms, neighborlist, atomic_weights, neighbor_indices, Seq)

where:
- center_atoms: Npairs × 3 matrix of central atom positions.
- neighborlist: Npairs × 3 matrix of displacement vectors (neighbour - center).
- atomic_weights: length Npairs vector of atomic numbers (signed if weight_on).
- neighbor_indices: Npairs × 2 matrix of [i, j] indices.
- Seq: Nseq × 2 matrix, for each i contains [i, i] and [i, j] for all neighbours j (sorted).

All positions are in Ångström (Float64).
"""
function build_neighbor_data(system, rcut::Real; weight_on::Bool=false)
    # 1. Extract atomic information
    atoms = collect(system)
    n_atoms = length(system)
    symbols = [String(atomic_symbol(atom)) for atom in atoms]
    atomic_nums = [ATOMIC_NUMBER[sym] for sym in symbols]
    # Convert positions to Float64 in Å
    pos_float = [ustrip.(u"Å", position(atom)) for atom in atoms]
    positions = reduce(vcat, [reshape(p, 1, 3) for p in pos_float])

    # 2. Get the bounding box (unit cell) and periodic flags
    cell = bounding_box(system)
    if cell === nothing
        error("System has no bounding box. Please provide a periodic system.")
    end
    cell_float = ustrip.(u"Å", cell)   # 3×3 matrix
    # Default to fully periodic if PBC flags are not stored
    pbc = haspbc(system) ? [pbc(system)...] : [true, true, true]

    # 3. Build neighbour list using NeighbourLists.jl
    nlist = neighbour_list(positions, rcut, cell_float, pbc)

    # 4. Iterate over all atoms and collect pair data
    center_atoms_list = Vector{Float64}[]
    neighborlist_list = Vector{Float64}[]
    atomic_weights_list = Int[]
    neighbor_indices_list = Vector{Int}[]

    for i in 1:n_atoms
        # Get neighbours: indices, displacement vectors (R), and images (not needed here)
        j_indices, R_vectors, _ = neighbours(nlist, i)

        center_pos = positions[i, :]
        for (j, R) in zip(j_indices, R_vectors)
            # R is the displacement from atom i to atom j (minimum image)
            push!(center_atoms_list, center_pos)
            push!(neighborlist_list, R)

            # Compute atomic weight (same logic as Python)
            if weight_on && atomic_nums[i] != atomic_nums[j]
                w = -atomic_nums[j]
            else
                w = atomic_nums[j]
            end
            push!(atomic_weights_list, w)
            push!(neighbor_indices_list, [i, j])
        end
    end

    # Convert lists to matrices (each row is one pair)
    center_atoms = reduce(vcat, [reshape(v, 1, 3) for v in center_atoms_list])
    neighborlist = reduce(vcat, [reshape(v, 1, 3) for v in neighborlist_list])
    neighbor_indices = reduce(vcat, [reshape(v, 1, 2) for v in neighbor_indices_list])
    atomic_weights = Int.(atomic_weights_list)

    # 5. Build Seq: for each atom i, include i itself and all unique neighbours (sorted)
    Seq_list = Vector{Int}[]
    for i in 1:n_atoms
        rows = findall(neighbor_indices[:, 1] .== i)
        unique_j = unique(neighbor_indices[rows, 2])
        push!(unique_j, i)          # self-pair
        sort!(unique_j)
        for j in unique_j
            push!(Seq_list, [i, j])
        end
    end
    Seq = reduce(vcat, [reshape(v, 1, 2) for v in Seq_list])

    return center_atoms, neighborlist, atomic_weights, neighbor_indices, Seq
end

# ============================
# Example usage
# ============================
filename = "DATASET_N2H4_v2/N2H4_2mol_1data.xyz"
system = read_first_frame(filename)

# Print atoms and lattice (matching the Python output)
for atom in system
    pos = ustrip.(u"Å", position(atom))
    println("$(atomic_symbol(atom)) $(pos[1]) $(pos[2]) $(pos[3])")
end
println("Lattice vectors: ")
cell = bounding_box(system)
if cell !== nothing
    display(ustrip.(u"Å", cell))
else
    println("No cell defined")
end

rcut = 3.5
weight_on = false   # default

center_atoms, neighborlist, atomic_weights, neighbor_indices, Seq =
    build_neighbor_data(system, rcut; weight_on=weight_on)

# Check sizes (optional)
@show size(center_atoms)
@show size(neighborlist)
@show size(atomic_weights)
@show size(neighbor_indices)
@show size(Seq)