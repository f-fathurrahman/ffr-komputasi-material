mutable struct Atoms
    Natoms::Int64
    positions::Array{Float64,2}
    energy::Float64
    forces::Array{Float64,2}
end

function Atoms(Natoms::Int64)
    positions = zeros(Float64,3,Natoms)
    energy = 0.0
    forces = zeros(Float64,3,Natoms)
    return Atoms(Natoms, positions, energy, forces)
end

import Base: length
function length(atoms::Atoms)
    return size(atoms.positions,2)
end