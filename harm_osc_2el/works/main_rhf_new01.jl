using Infiltrator
using LinearAlgebra: Diagonal, kron

using MyModule

function main()
    n = 2
    l = 20

    ω = 0.25
    basis = HarmonicOscillatorBasis(l, ω) |> SpinBasis
    V = HarmonicOscillatorCoulomb(ω, shielding=0.25)

    xgrid = range(-10, stop=10, length=2001) |> collect

    @infiltrate
    # Debug SpatialSystem
    #system = init_system(n, basis, xgrid, V)

end

main()