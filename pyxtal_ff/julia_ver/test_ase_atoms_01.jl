using LinearAlgebra
using StaticArrays

mutable struct Atoms
    symbols::Vector{String}
    positions::Vector{SVector{3, Float64}}
    cell::SMatrix{3, 3, Float64, 9}  # SMatrix is great here too since it's always 3x3
    pbc::SVector{3, Bool}
    tags::Vector{Int}
    masses::Vector{Float64}
    info::Dict{String, Any}

    function Atoms(symbols, positions, cell, pbc, tags, masses, info)
        n = length(symbols)
        length(positions) == n || throw(DimensionMismatch("positions length must match symbols"))
        new(symbols, positions, cell, pbc, tags, masses, info)
    end
end


# Ergonomic outer constructor
function Atoms(;
    symbols::AbstractVector{String},
    positions::AbstractVector{<:Any}, # Accepts Vector of Vectors, Tuples, or SVectors
    cell::AbstractMatrix{<:Real} = zeros(3, 3),
    pbc::Union{AbstractVector{Bool}, Bool} = [false, false, false],
    tags::AbstractVector{Int} = zeros(Int, length(symbols)),
    masses::AbstractVector{<:Real} = zeros(Float64, length(symbols)),
    info::Dict{String, Any} = Dict{String, Any}()
)
    # Convert inputs to static types for performance
    pos_svec = [SVector{3, Float64}(p) for p in positions]
    cell_smat = SMatrix{3, 3, Float64}(cell)
    pbc_svec = pbc isa Bool ? SVector{3, Bool}(pbc, pbc, pbc) : SVector{3, Bool}(pbc...)
    
    return Atoms(collect(symbols), pos_svec, cell_smat, pbc_svec, collect(tags), collect(masses), info)
end

# Basic properties
Base.length(atoms::Atoms) = length(atoms.symbols)
get_positions(atoms::Atoms) = atoms.positions
get_cell(atoms::Atoms) = atoms.cell
get_pbc(atoms::Atoms) = atoms.pbc
get_tags(atoms::Atoms) = atoms.tags

# Geometry and Box operations
function get_volume(atoms::Atoms)
    return abs(det(atoms.cell))
end

# Slicing/Indexing to extract subsets (e.g., getting just the adsorbate)
function Base.getindex(atoms::Atoms, idx)
    return Atoms(
        symbols = atoms.symbols[idx],
        positions = atoms.positions[idx, :],
        cell = atoms.cell,
        pbc = atoms.pbc,
        tags = atoms.tags[idx],
        masses = atoms.masses[idx],
        info = copy(atoms.info)
    )
end

# Display override for clean REPL output
function Base.show(io::IO, atoms::Atoms)
    formula = join(unique(atoms.symbols)) # Simplified formula logic
    print(io, "Atoms(symbols=\"$(formula)\", n_atoms=$(length(atoms)), pbc=$(atoms.pbc))")
end

# Translation is incredibly clean
function translate!(atoms::Atoms, translation::AbstractVector{<:Real})
    t_vec = SVector{3, Float64}(translation)
    atoms.positions .+= Ref(t_vec) # Ref prevents broadcasting over the 3 elements of t_vec
    return atoms
end

# Cell transformation using broadcasting
function set_cell!(atoms::Atoms, new_cell::AbstractMatrix{<:Real}; scale_atoms::Bool=false)
    new_cell_s = SMatrix{3, 3, Float64}(new_cell)
    
    if scale_atoms
        inv_cell = inv(atoms.cell)
        # Apply transformation to every SVector in the array
        atoms.positions = map(atoms.positions) do pos
            fractional = inv_cell * pos
            new_cell_s * fractional
        end
    end
    atoms.cell = new_cell_s
end


# Concatenation naturally uses standard vcat for the arrays
function Base.vcat(a1::Atoms, a2::Atoms)
    return Atoms(
        symbols = vcat(a1.symbols, a2.symbols),
        positions = vcat(a1.positions, a2.positions), # standard vcat works for Vector{SVector}
        cell = a1.cell,
        pbc = a1.pbc,
        tags = vcat(a1.tags, a2.tags),
        masses = vcat(a1.masses, a2.masses),
        info = copy(a1.info)
    )
end


function test_main()

    # 1. Create a simple 2-atom Pt slab
    slab = Atoms(
        symbols=["Pt", "Pt"],
        positions=[
            [0.0, 0.0, 0.0], 
            [1.5, 1.5, 0.0]
        ],
        cell=[3.0  0.0  0.0; 
            0.0  3.0  0.0; 
            0.0  0.0 15.0], 
        pbc=[true, true, false],                       
        tags=[1, 1]                                    
    )

    # 2. Create the CO adsorbate
    co_molecule = Atoms(
        symbols=["C", "O"],
        positions=[
            [0.0, 0.0, 0.0], 
            [0.0, 0.0, 1.128]
        ],
        pbc=[false, false, false],
        tags=[2, 2]                                    
    )

    # Move the CO molecule 2.0 Å above the surface
    translate!(co_molecule, [0.0, 0.0, 2.0])

    # 3. Combine them using standard vcat
    system = vcat(slab, co_molecule)

    println(system)
    # Atoms(symbols="PtCO", n_atoms=4, pbc=Bool[1, 1, 0])

    # Grab the positions array
    pos = get_positions(system)

    # Display the whole array of SVectors
    display(pos)
    # 4-element Vector{SVector{3, Float64}}:
    #  [0.0, 0.0, 0.0]
    #  [1.5, 1.5, 0.0]
    #  [0.0, 0.0, 2.0]
    #  [0.0, 0.0, 3.128]

    # Calculate the distance between Pt(2) and C(3) with zero heap allocations
    distance = norm(pos[2] - pos[3])
    println("Pt-C distance: ", round(distance, digits=3), " Å")
    # Pt-C distance: 2.915 Å

    return
end

