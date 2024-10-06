import pickle
import numpy as np
from draw_strat_sample import *

np.random.seed(1234)

with open("ALL_ATOMS.pkl", "rb") as f:
    all_atoms = pickle.load(f)

Ndata = len(all_atoms)
Natoms = len(all_atoms[0])

energies = np.zeros((Ndata))
for i in range(Ndata):
    energies[i] = all_atoms[i].calc.get_potential_energy()

Ntrain = 2000
idx_train = draw_strat_sample(energies, Ntrain)
energies_train = energies[idx_train]

boxes_train = np.zeros((Ntrain, 3, 3))
forces_train = np.zeros((Ntrain, Natoms, 3))
coords_train = np.zeros((Ntrain, Natoms, 3))
for i,itrain in enumerate(idx_train):
    coords_train[i,:,:] = all_atoms[itrain].get_positions()
    forces_train[i,:,:] = all_atoms[itrain].calc.get_forces()
    boxes_train[i,:,:] = all_atoms[itrain].get_cell()[:]
    all_atoms[itrain].write("TRAIN.xyz", append=True)

np.save("box.npy", boxes_train)
np.save("coord.npy", coords_train)
np.save("force.npy", forces_train)
np.save("energy.npy", energies_train)
np.save("atomic_number.npy", all_atoms[0].get_atomic_numbers())


