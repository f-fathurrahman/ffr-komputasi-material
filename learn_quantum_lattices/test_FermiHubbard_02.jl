Pkg.activate("quantum_lattices", shared=true)

push!(LOAD_PATH, ".")

using MyQuantumLattices
using MyExactDiagonalization
using LinearAlgebra: eigen

lattice = Lattice([0.0], [1.0])

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

bases = BinaryBases(1:2, 1)

t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 8.0)

#=
# define the Hilbert space (single-orbital spin-1/2 complex fermion)
hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(unitcell))

# define the binary bases of the a half-filled system on the above cluster
bases = BinaryBases(1:12, 6) ⊗ BinaryBases(13:24, 6)

# define the terms, i.e. the nearest-neighbor hopping and the Hubbard interaction
t = Hopping(:t, -1.0, 1)
#display(t)

U = Hubbard(:U, 8.0)
#display(U)

# define the exact diagonalization algorithm for the Fermi Hubbard model
ed = ED(lattice, hilbert, (t, U), TargetSpace(bases))

# find the ground state and its energy
eigensystem = eigen(matrix(ed); nev=1)

# Ground state energy should be -4.913259209075605
println(eigensystem.values)
=#
