Pkg.activate("quantum_lattices", shared=true)

push!(LOAD_PATH, ".")

using MyQuantumLattices
using MyExactDiagonalization
using LinearAlgebra: eigen

@inline function debug_matrix(
    ed::ED,
    braket::NTuple{2, Sector}=first(ed.Hₘ.transformation.brakets);
    kwargs...
)
    return expand(SectorFilter(braket)(ed.Hₘ))[braket]
end



lattice = Lattice([0.0], [1.0])

hilbert = Hilbert(site=>Fock{:f}(1, 2) for site=1:length(lattice))

#bases = BinaryBases(1:2, 1) ⊗ BinaryBases(3:4, 1)
bases = BinaryBases(3:4, 1) ⊗ BinaryBases(1:2, 1)
display(bases); println()

t = Hopping(:t, -1.0, 1)
U = Hubbard(:U, 8.0)

ed = ED(lattice, hilbert, (t, U), TargetSpace(bases))

Ham1 = matrix(ed)
#display(Ham1.matrix)

println("Pass here")