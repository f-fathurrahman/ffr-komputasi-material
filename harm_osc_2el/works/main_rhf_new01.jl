using Infiltrator
using LinearAlgebra: Diagonal, kron

include("ffr/HarmonicOscillatorBasis.jl")
include("ffr/interactions.jl")


function calc_onebody_integrals(ho::HarmonicOscillatorBasis, xgrid)
    l = ho.l
    ω = ho.ω
    return Diagonal([(n + 0.5) * ω for n in 0:l-1])
end
# xgrid is not used?


function calc_onebody_integrals(basis::SpinBasis, xgrid)
    h = calc_onebody_integrals(basis.base, xgrid)
    h = kron(h, [1 0; 0 1])
    return h
end

function calc_twobody_integrals(spfs, grid, V::Interaction)
    inner = _twobody_inner_ints(spfs, grid, V)
    u = _twobody_outer_int(spfs, grid, inner)
    return u
end

function _twodbody_inner_ints(spfs, grid, V::Interaction)
    l = length(spfs)
    inner_int = zeros(typeof(spfs[1][1]), l, l, length(spfs[1]))
    
    #interactions = [similar(grid) for i in 1:Threads.nthreads()]
    #fs = [similar(spfs[1]) for i in 1:Threads.nthreads()]
    
    for xi in eachindex(grid)
        x1 = grid[xi]
        #f_vals = fs[Threads.threadid()] # Pre-allocated vector for this thread
        #interaction = interactions[Threads.threadid()] # ^^^
        evaluate_on_grid!(interaction, x1, grid, V)
        for κ in 1:l
            for λ in κ:l
                f_vals .= conj.(spfs[κ]) .* interaction .* spfs[λ]
                res = trapz(f_vals, grid)
                inner_int[κ, λ, xi] = res
                inner_int[λ, κ, xi] = res'
            end
        end
    end
    return inner_int
end



function main()
    n = 2
    l = 20

    ω = 0.25
    basis = HarmonicOscillatorBasis(l, ω) |> SpinBasis
    V = HarmonicOscillatorCoulomb(ω, shielding=0.25)

    xgrid = range(-10, stop=10, length=2001)
    
    # Debug SpatialSystem

    # The basis functions evaluated on the grid
    spfs = evaluate_on_grid(basis, xgrid)
    h = calc_onebody_integrals(basis, xgrid) # One body integrals

    @infiltrate
end

main()