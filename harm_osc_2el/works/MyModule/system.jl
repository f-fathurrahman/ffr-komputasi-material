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
# XXX Pass SpatialBasis ?
# SpatialBasis is an abstract type
function init_system(n, basis::SpatialBasis, xgrid, V::Interaction)
    Nbasis = basis.l
    spfs = evaluate_on_grid(basis, xgrid) # including spin
    h = calc_onebody_integrals(basis, xgrid) # One body integrals
    u = calc_twobody_integrals(basis, xgrid, V) # Two body integrals
    u .= u .- permutedims(u, [1, 2, 4, 3]) # Anti-symmetrizing u
    transform = LinearAlgebra.I(Nbasis)
    return SpatialSystem{typeof(basis)}(
        n, Nbasis, h, u, spfs, xgrid, basis, transform, V
    )
end

