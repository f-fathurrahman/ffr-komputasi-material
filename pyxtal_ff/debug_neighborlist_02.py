import ase.io
from ase.neighborlist import NeighborList
import numpy as np

atoms_list = ase.io.read("DATASET_OTHERS/TiAl_gabung.xyz@:")
#atoms_list = ase.io.read("DATASET_N2H4_v2/N2H4_2mol_1data.xyz@:")
#atoms_list = ase.io.read("DATASET_N2H4_v1/TEMP_ATOMS_TRAIN.xyz@:")
atoms = atoms_list[0]
for a in atoms:
    print(f"{a.symbol} {a.position[0]} {a.position[1]} {a.position[2]}")

print("Lattice vectors: ")
print(atoms.get_cell())

rcut = 3.5

# default
weight_on = False

atom_ids = range(len(atoms))

cutoffs = [rcut/2]*len(atoms)

# self.primitive is assumed to be false

nl = NeighborList(cutoffs, self_interaction=False, bothways=True, skin=0.0)
nl.update(atoms)

center_atoms = []
neighbors = []
neighbor_indices = []
atomic_weights = []
temp_indices = []

# Loop over atoms
for i in atom_ids:
    # get center atom position vector
    center_atom = atoms.positions[i]
    print(f"\nCenter atom: ia={i+1} {atoms[i].symbol}, {atoms[i].position}")
    # get indices and cell offsets for each neighbor
    indices, offsets = nl.get_neighbors(i)
    # Loop over neighbors
    for j, offset in zip(indices, offsets):
        print(f"idx={j} offset={offset}")
        # Get positions of the neighbor
        pos = atoms.positions[j] + np.dot(offset,atoms.get_cell())
        print("pos = ", pos)
        dist_vec = pos - center_atom
        print("dist_vec = ", dist_vec, " norm = ", np.linalg.norm(dist_vec))
