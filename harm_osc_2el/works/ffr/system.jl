abstract type System end

struct SpatialSystem{T} <: System
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    #
    h::Array{Float64, 2} # one-body integrals
    u::Array{Float64, 4} # two-body integrals
    spfs::Vector{Vector{Float64}} # the basis functions evaluated on the grid
    #
    grid::Vector{Float64}
    basis::T
    transform::Matrix{Float64}
    V::Interaction
end

# A constructor for SpatialSystem
function SpatialSystem(n, basis::SpatialBasis, grid, V::Interaction)
    l = basis.l
    spfs = spatial(basis, grid) # The basis functions evaluated on the grid
    h = onebody(basis, grid) # One body integrals
    u = twobody(basis, grid, V) # Two body integrals
    u .= u .- permutedims(u, [1, 2, 4, 3]) # Anti-symmetrizing u
    transform = LinearAlgebra.I(l)
    return SpatialSystem{typeof(basis)}(n, l, h, u, spfs, grid, basis, transform, V)
end

