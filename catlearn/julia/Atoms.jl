mutable struct Atoms
    Natoms::Int64
    Nspecies::Int64
    positions::Array{Float64,2}
    atm2species::Array{Int64,1}
    atsymbs::Array{String,1}         # for each atom
    SpeciesSymbols::Array{String,1}  # unique symbols
    LatVecs::Array{Float64,2}
    Zvals::Array{Float64,1}   # unique
    masses::Array{Float64,1}
    energy::Float64
    forces::Array{Float64,2}
    pbc::Tuple{Bool,Bool,Bool}
end

function Atoms(Natoms::Int64)
    Nspecies = 1
    positions = zeros(Float64,3,Natoms)
    atm2species = [1]
    atsymbs = ["X"]
    SpeciesSymbols = ["X"]
    LatVecs = diagm([1.0, 1.0, 1.0])
    Zvals = [1.0]
    masses = [1.0]
    energy = 0.0
    forces = zeros(Float64,3,Natoms)
    pbc = (false,false,false)
    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols,
        LatVecs, Zvals, masses, energy, forces, pbc)
end

