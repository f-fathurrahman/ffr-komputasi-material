using HarmOsc2El

n = 2 # not required for basis
l = 20 # This is the number of basis functions
ω = 0.25 # harmonic potential strength from which the basis functions are consctructed
basis = HOBasis(l, ω)
