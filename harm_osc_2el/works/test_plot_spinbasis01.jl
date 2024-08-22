using HarmOsc2El

import Plots
import PlotThemes
Plots.theme(:dark)

l = 20
ω = 0.25
basis = SpinBasis(HOBasis(l, ω))

# Potetial
V = HOCoulomb(ω, shielding = 0.25)

# Grid points
n = 2
Npoints = 2001
xgrid = range(-10, stop=10, length=2001) |> collect

χ = HarmOsc2El.evaluate(basis, xgrid)

Plots.plot(xgrid, χ[1]) # χ[1] == χ[2], spin degenerate
Plots.plot(xgrid, χ[3])
Plots.plot(xgrid, χ[5])