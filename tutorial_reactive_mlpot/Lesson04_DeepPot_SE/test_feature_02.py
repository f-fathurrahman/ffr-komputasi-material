import math
from typing import Sequence, Tuple

import numpy as np

import torch
from torch import Tensor
from torch import nn

from load_data import load_data
from models import *
from local_environment import local_environment
from  feature import Feature

coord, atom_types, grad = load_data()

descriptor = Feature(4, neuron=[25, 50], axis_neuron=4)

"""
# debug forward: pass coord and atom_types
coords = coord #[2:3]
atom_types = atom_types #[2:3]

num_batches, num_channels, _ = coords.size()
loc_env_r, loc_env_a = local_environment(coords)

neighbor_types = atom_types.repeat(num_channels, 1)
mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device)
neighbor_types = torch.masked_select(neighbor_types, mask).view(num_channels, -1)
indices = ((atom_types * descriptor.n_types).unsqueeze(-1) + neighbor_types).view(-1)

output, _ = descriptor.local_embedding((loc_env_r.view(num_batches, -1, 1), indices))
output = output.view(num_batches, num_channels, num_channels - 1, -1)

output = torch.transpose(output, 2, 3) @ (loc_env_a @ (torch.transpose(loc_env_a, 2, 3) @ output[..., :descriptor.axis_neuron]))
output = output.view(num_batches, num_channels, -1)
"""