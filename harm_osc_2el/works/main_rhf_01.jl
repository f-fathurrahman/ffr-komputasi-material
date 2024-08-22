using HarmOsc2El

n = 2
l = 20

ω = 0.25
basis = SpinBasis(HOBasis(l, ω))

# Potetial
V = HOCoulomb(ω, shielding = 0.25)

# Grid points
grid = [x for x in range(-10, stop = 10, length = 2001)]

#=
Construct system with the following inputs:
- number of points: n (redundant? can be obtained from grid?)
- basis: spin basis
- grid points: spatial points array
- potential: V
=#
system = System(n, basis, grid, V)

# Construct various arrays needed in RHF calculations
rhf = RHF(system)
compute_ground_state!(rhf);

# XXX Why need this?
rhf_system = System(rhf)
println("Reference energy: $(reference_energy(rhf_system)), $(t) s")

