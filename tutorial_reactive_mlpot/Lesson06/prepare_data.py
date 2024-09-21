import numpy as np
import torch

"""
unit           conversion factor
1 Hartree       627.509 kcal mol-1
1 Hartree	    27.2116 eV
1 kcal mol-1    4.184 kJ mol-1
1 eV            23.06035 kcal mol-1
"""

def prepare_data():

    qm_coord = np.array(np.load("../DATASET/Claisen_Rearrangement/qm_coord.npy"), dtype="float32")
    atom_types = np.loadtxt("../DATASET/Claisen_Rearrangement/qm_elem.txt", dtype=int)
    elems = np.unique(atom_types).tolist()
    atom_types = np.array([[elems.index(i) for i in atom_types]])
    atom_types = atom_types.repeat(len(qm_coord), axis=0)

    energy = np.array(
        (np.load("../DATASET/Claisen_Rearrangement/energy.npy") -
        np.load("../DATASET/Claisen_Rearrangement/energy_sqm.npy")
        ) * 27.2114 * 23.061, dtype="float32"
    )
    energy = energy - energy.mean()
    qm_gradient = np.array(
        (np.load("../DATASET/Claisen_Rearrangement/qm_grad.npy") -
        np.load("../DATASET/Claisen_Rearrangement/qm_grad_sqm.npy")
        ) * 27.2114 * 23.061 / 0.529177249, dtype="float32"
    )

    device = "cpu"

    qm_coord = torch.from_numpy(qm_coord).to(device)
    atom_types = torch.from_numpy(atom_types).to(device)
    energy = torch.from_numpy(energy).to(device)
    qm_gradient = torch.from_numpy(qm_gradient).to(device)

    return qm_coord, atom_types, energy, qm_gradient
