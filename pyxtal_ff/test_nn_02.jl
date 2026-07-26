using AtomsIO, Unitful, NeighbourLists

function read_first_frame(filename::String)
    # load_system returns an AbstractSystem with the first frame
    system = load_system(filename)
    return system
end

function debug_main()
    filename = "DATASET_N2H4_v2/N2H4_2mol_1data.xyz"
    system = load_system(filename, 1) # first frame

    # Print atoms and lattice (matching the Python output)
    #for atom in system
    #    pos = ustrip.(u"Å", position(atom))
    #    println("$(atomic_symbol(atom)) $(pos[1]) $(pos[2]) $(pos[3])")
    #end
    #cell_vectors = system.system_data.cell_vectors

    atoms = collect(system)
    Natoms = length(system)
    symbols = [atomic_symbol(atom) for atom in atoms]
    atomic_nums = [atomic_number(atom) for atom in atoms]

    #=
    cell_vectors = [zeros(3) for _ in 1:3]
    for i in 1:3, j in 1:3
        cell_vectors[i][j] = ustrip(system.system_data.cell_vectors[i][j])
    end
    =#
    cell_vectors = zeros(Float64, 3, 3)
    for i in 1:3, j in 1:3
        # XXX Is this the actual indexing?
        cell_vectors[i,j] = ustrip(system.system_data.cell_vectors[i][j])
    end

    positions = Vector{SVector{3, Float64}}(undef,Natoms)
    for ia in 1:Natoms
        r = position(atoms[ia])
        x = ustrip(r[1])
        y = ustrip(r[2])
        z = ustrip(r[3])
        positions[ia] = @SVector [x,y,z]
    end

    # nlist = neighbour_list(system, 3.5u"Å") # directly using AtomsBase
    # 3. Build neighbour list using NeighbourLists.jl
    pbc = [true, true, true]
    nlist = neighbour_list(positions, 3.5, cell_vectors, pbc)

    center_atoms_list = Vector{Float64}[]
    neighborlist_list = Vector{Float64}[]
    atomic_weights_list = Int[]
    neighbor_indices_list = Vector{Int}[]

    for i in 1:Natoms
        # Get neighbours: indices, displacement vectors (R), and images (not needed here)
        j_indices, R_vectors, _ = neighbours(nlist, i)
        println("j_indices = ", j_indices)
    end

    @exfiltrate

    return
end

