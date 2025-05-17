import ase.io
from ase.neighborlist import NeighborList
import numpy as np

#atoms_list = ase.io.read("DATASET_OTHERS/TiAl_gabung.xyz@:")
#atoms_list = ase.io.read("DATASET_N2H4_v2/N2H4_2mol_1data.xyz@:")
atoms_list = ase.io.read("DATASET_N2H4_v1/TEMP_ATOMS_TRAIN.xyz@:")
atoms = atoms_list[0]
for a in atoms:
    print(f"{a.symbol} {a.position[0]} {a.position[1]} {a.position[2]}")

print("Lattice vectors: ")
print(atoms.get_cell())


# Should be invariant with w.r.t translations
#pos_shifted = atoms.positions.copy()
#pos_shifted[:,2] = atoms.positions[:,2] + 3.0
#atoms.set_positions(pos_shifted)

lmax = 4
nmax = 3
rcut = 3.5
alpha = 2.0

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
    temp_indices.append(indices)
    # Loop over neighbors
    for j, offset in zip(indices, offsets):
        print(f"idx={j} offset={offset}")
        # wrap positions to unit cell
        pos = atoms.positions[j] + np.dot(offset,atoms.get_cell()) - center_atom
        print(f"distance vector within unit cell = {pos}")
        center_atoms.append(center_atom) # current center atom
        neighbors.append(pos)
        if weight_on and atoms[j].number != atoms[i].number:
            factor = -1
        else:
            factor = 1
        atomic_weights.append(factor*atoms[j].number)
        neighbor_indices.append([i,j]) # the pair

# convert neighbor_indices to np.array
neighbor_indices = np.array(neighbor_indices, dtype=np.int64)

Seq = []
for i in atom_ids:
    ineighs = neighbor_indices[:,0] == i
    print(f"\nneighbor_indices for atom {i} = {neighbor_indices[ineighs]}")
    unique_atoms = np.unique(neighbor_indices[ineighs]) # unique atom indices ?
    print("unique atoms = ", unique_atoms)
    #
    # Handle the case i is not in unique neighbor atoms
    if i not in unique_atoms:
        # XXX append and sort?
        at = list(unique_atoms)
        at.append(i)
        at.sort()
        unique_atoms = np.array(at)
    print("After if unique_atoms = ", unique_atoms)
    for j in unique_atoms:
        Seq.append([i,j])
    #print("Current Seq = ", Seq)

# These are set internal variables
Seq = np.array(Seq, dtype=np.int64)
center_atoms = np.array(center_atoms, dtype=np.float64)
neighborlist = np.array(neighbors, dtype=np.float64)
seq = Seq
atomic_weights = np.array(atomic_weights, dtype=np.int64)
neighbor_indices = neighbor_indices
