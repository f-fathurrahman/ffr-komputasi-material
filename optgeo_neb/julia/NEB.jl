mutable struct NEBCalculator
    images::Vector{Atoms}
    k::Vector{Float64} # spring constant    
    Nimages::Int64
    Natoms::Int64
    pbc::Tuple{Bool,Bool,Bool}
    energies::Vector{Float64}
    neb_forces::Array{Float64,3}
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
    neb_forces = zeros(Float64,3,Natoms,Nimages-2)
    #
    k = 0.1*ones(Nimages-1)
    #
    interpolate!(images)
    return NEBCalculator(images, k, Nimages, Natoms, pbc, energies, neb_forces)
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

function setup_initial_final!(calc, neb::NEBCalculator)
    calc_energy_forces!(calc, neb.images[1])
    calc_energy_forces!(calc, neb.images[end])
    return
end

function compute!(calc, neb::NEBCalculator)
    
    # forces = zeros(3,Natoms,Nimages-2)
    
    forces = neb.neb_forces
    Nimages = neb.Nimages
    images = neb.images
    energies = neb.energies
    k = neb.k

    for i in 2:Nimages-1
        calc_energy_forces!(calc, images[i])
        # copy results
        energies[i] = images[i].energy
        forces[:,:,i-1] = images[i].forces[:,:]
        println("energy = ", energies[i])
        println("forces = ", forces[:,:,i-1])
    end

    t1 = images[2].positions - images[1].positions
    nt1 = norm(t1)

    for i in 2:Nimages-1
        t2 = images[i+1].positions - images[i].positions
        nt2 = norm(t2)
        tangent = t1 + t2
        tt = dot(tangent, tangent)
        #
        @views f = forces[:,:,i-1]
        ft = dot(f, tangent)
        f = f - ft / tt * tangent
        f = f - dot(t1*k[i-1] - t2*k[i], tangent) / tt * tangent
        #
        t1 = t2
        nt1 = nt2
        println("forces = ", forces[:,:,i-1])
    end
    return
end