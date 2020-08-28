mutable struct NEBCalculator
    images::Vector{Atoms}
    Nimages::Int64
    Natoms::Int64
    pbc::Tuple{Bool,Bool,Bool}
    energies::Vector{Float64}
    real_forces::Vector{Matrix{Float64}}
end

function NEBCalculator(images::Vector{Atoms})
    Nimages = length(images)
    # Number of atoms
    Natoms = images[1].Natoms
    for i in 2:Nimages
        @assert Natoms == images[i].Natoms
    end
    # PBC or not
    pbc = images[1].pbc
    for i in 2:Nimages
        @assert pbc == images[i].pbc
    end
    #
    energies = zeros(Float64,Nimages)
    real_forces = Vector{Matrix{Float64}}(undef,Nimages)
    for i in 1:Nimages
        real_forces[i] = zeros(Float64,3,Natoms)
    end
    #
    interpolate!(images)
    return NEBCalculator(images, Nimages, Natoms, pbc, energies, real_forces)
end

# Linear interpolation
function interpolate!( images::Vector{Atoms} )
    Nimages = length(images)
    pos1 = images[1].positions
    pos2 = images[Nimages].positions
    d = (pos2 - pos1)/(Nimages-1)
    #
    # FIXME: Consider minimum image convention
    #
    for i in 2:Nimages
        images[i].positions[:] = pos1[:] + (i-1)*d[:]
    end
    return
end
