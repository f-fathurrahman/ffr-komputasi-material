from ase.build import bulk
from my_neighborlist import MyNeighborList

a_lattice = 4.05  # Angstrom lattice spacing
atoms = bulk("Al", "fcc", a=a_lattice)

rc_list = 2.4
nl = MyNeighborList([0.5 * rc_list] * len(atoms), self_interaction=False)
nl.update(atoms)
print("Number of neighbors: ", nl.nneighbors)

idx_atom = 0
print("Position of atom: ", atoms.positions[idx_atom])
indices, offsets = nl.get_neighbors(idx_atom)
for i, offset in zip(indices, offsets):
    print(atoms.positions[i] + offset @ atoms.get_cell())

