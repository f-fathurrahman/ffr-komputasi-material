from typing import Sequence

import numpy as np

import torch
from torch import Tensor

def load_data():
    ds = np.DataSource(None)
    coord = np.array(np.load(ds.open("../DATASET/DeepPot/input_coord.npy", "rb")), dtype="float32")
    atom_types = np.loadtxt(ds.open("../DATASET/DeepPot/type.raw", "r"), dtype=int)
    elems = np.unique(atom_types).tolist()
    atom_types = np.array([[elems.index(i) for i in atom_types]])
    atom_types = atom_types.repeat(len(coord), axis=0)
    grad = np.array(np.load(ds.open("../DATASET/DeepPot/input_grad.npy", "rb")), dtype="float32")

    # Now, convert them to Tensor
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    coord = torch.from_numpy(coord).to(device)
    atom_types = torch.from_numpy(atom_types).to(device)
    grad = torch.from_numpy(grad).to(device)
    return coord, atom_types, grad


coord, atom_types, grad = load_data()


coords = coord[0:1] # a single molecule (one data point)


num_batches, num_channels, _ = coords.size()
# num_channels: no. of atoms (?)

# computes pairwise distance vectors
rij_ = coords[:, :, None] - coords[:, None]
dij_ = torch.norm(rij_, dim=3)

mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device) # remove self-interaction
# inverse of unit matrix of type bool
# all diagonal elements are False
# other elements are True

rij = torch.masked_select(rij_, mask.unsqueeze(2)).view(num_batches, num_channels, num_channels - 1, 3)
dij = torch.masked_select(dij_, mask).view(num_batches, num_channels, num_channels - 1)

dij_inv = 1 / dij
dij2_inv = dij_inv * dij_inv

loc_env_r = dij_inv  # distances
loc_env_a = rij * dij2_inv.unsqueeze(3)  # vectors

