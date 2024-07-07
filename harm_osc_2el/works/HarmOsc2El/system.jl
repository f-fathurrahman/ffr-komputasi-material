abstract type System end

include("BasisFunctions.jl")

include("interactions.jl")
include("spatialintegrals.jl")
include("spinintegrals.jl")

struct SpatialSystem{T} <: System
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    
    h::Array{Float64, 2} # one-body integrals
    u::Array{Float64, 4} # two-body integrals
    spfs::Vector{Vector{Float64}} # the basis functions evaluated on the grid
    
    grid::Vector{Float64}
    basis::T
    transform::Matrix{Float64}
    V::Interaction
end

function System(n, basis::SpatialBasis, grid, V::Interaction)
    l = basis.l
    spfs = spatial(basis, grid) # The basis functions evaluated on the grid
    h = onebody(basis, grid) # One body integrals
    u = twobody(basis, grid, V) # Two body integrals
    u .= u .- permutedims(u, [1, 2, 4, 3]) # Anti-symmetrizing u

    transform = LinearAlgebra.I(l)
    return SpatialSystem{typeof(basis)}(n, l, h, u, spfs, grid, basis, transform, V)
end

function reference_energy(system)
    (; h, n, u) = system
    E = 0.0
    @inbounds for i in 1:n
        E += h[i, i]
        for j in 1:n
            E += 0.5 * u[i, j, i, j]
        end
    end
    return E
end

function sp_energies(system)
    (; l) = system
    
    f = fock_matrix(system)
    ϵ = zeros(l)
    @inbounds for q in 1:l
        ϵ[q] = f[q, q]
    end
    return ϵ
end

function fock_matrix(system::System)
    (; n, l, h, u) = system


    P = Float64.(la.diagm( vcat( repeat([1], n) , repeat([0], l-n) ) ))
    
    F = zeros(l, l)
    F += h
    for c in 1:l
        for d in 1:l
            @inbounds P_dc = P[d, c]
            for a in 1:l
                for b in 1:l
                    @inbounds F[a, b] += P_dc * u[a, c, b, d]
                end
            end
        end
    end
    return F
end

struct PairingSystem <: System
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    
    h::Array{Float64, 2} # one-body integrals
    u::Array{Float64, 4} # two-body integrals
end

function System(n, basis::Pairing)
    (; l, states) = basis
    
    h = zeros((l, l))
    u = zeros((l, l, l, l))
    for p in 1:l
        for q in 1:l
            # One body integrals
            h[p, q] = Ĥ₀(states[p], states[q])

            for r in 1:l
                for s in 1:l
                    # Two body integrals
                    u[p, q, r, s] = V̂(states[p], states[q], states[r], states[s])
                end
            end
        end
    end

    return PairingSystem(n, l, h, u)
end

function sp_energies(system::PairingSystem)
    (; l, h, u, n) = system
    
    ϵ = zeros(l)
    @inbounds for q in 1:l
        ϵ[q] = h[q, q]
        for i in 1:n
            ϵ[q] += u[q, i, q, i]
        end
    end
    return ϵ
end

function particle_density(system)
    (;spfs, n) = system
    
    p_density = zero(spfs[1])
    for i in 1:n
        p_density .+= spfs[i].^2
    end
    return p_density
end

include("transform.jl")
;