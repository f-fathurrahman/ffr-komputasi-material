using Infiltrator
using LinearAlgebra: Diagonal, kron

using MyModule

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