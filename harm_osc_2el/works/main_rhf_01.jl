using HarmOsc2El

n = 2
l = 20

ω = 0.25
basis = SpinBasis(HOBasis(l, ω))
# Potetial (using the same ω as used in constructing the basis)
V = HOCoulomb(ω, shielding = 0.25)

# Grid points
xgrid = range(-10, stop=10, length=2001) |> collect

#=
Construct system with the following inputs:
- number of particles (electrons): n
- basis: spin basis
- grid points: spatial points array
- potential: V
=#
system = System(n, basis, xgrid, V)
# system of type SpatialSystem (?)

# Construct various arrays needed in RHF calculations
rhf = RHF(system)
t = @elapsed compute_ground_state!(rhf);

# XXX Why need this?
rhf_system = System(rhf)
println("Reference energy: $(reference_energy(rhf_system)), $(t) s")

