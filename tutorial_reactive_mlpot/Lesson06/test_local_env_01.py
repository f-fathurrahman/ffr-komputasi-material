from typing import Sequence
import numpy as np
import torch
from torch import Tensor
from prepare_data import *

qm_coord, atom_types, energy, qm_gradient = prepare_data()

Ndata = 2
coords = qm_coord[0:Ndata]
# coords have shape (Ndata,Natoms,3)
# Ndata will be referred to num_batches
# Natoms is num_channels

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
# rij contains the vectors last dimension is 3

dij = torch.masked_select(dij_, mask).view(num_batches, num_channels, num_channels - 1)
# dij is the magnitude of the vector

dij_inv = 1 / dij
dij2_inv = dij_inv * dij_inv

loc_env_r = dij_inv  # distances
loc_env_a = rij * dij2_inv.unsqueeze(3)  # vectors

