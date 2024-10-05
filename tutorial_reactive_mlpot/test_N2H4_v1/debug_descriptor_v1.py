import numpy as np
import torch
from my_deeppot import *
from prepare_data import *

torch.manual_seed(1234)

def debug_feature_forward(feature, coords, atom_types):
    num_batches, num_channels, _ = coords.size()
    loc_env_r, loc_env_a = local_environment(coords)

    n_types = feature.n_types
    local_embedding = feature.local_embedding
    axis_neuron = feature.axis_neuron

    neighbor_types = atom_types.repeat(num_channels, 1)
    mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device)
    neighbor_types = torch.masked_select(neighbor_types, mask).view(num_channels, -1)
    indices = ((atom_types * n_types).unsqueeze(-1) + neighbor_types).view(-1)

    output, _ = local_embedding((loc_env_r.view(num_batches, -1, 1), indices))
    output = output.view(num_batches, num_channels, num_channels - 1, -1)

    output = torch.transpose(output, 2, 3) @ (loc_env_a @ (torch.transpose(loc_env_a, 2, 3) @ output[..., :axis_neuron]))
    output = output.view(num_batches, num_channels, -1)

    return output


## XXX atom_types_all is not used in any of Feature, Fitting, or DeepPot instances

qm_coord, atom_types_all, energy, qm_gradient = prepare_data()

# coordinate, only one data point
coords = qm_coord[0:1,:,:]
atom_types = atom_types_all[0]

#coords.requires_grad_() # no need to track the gradient

#descriptor = Feature(3, neuron=[25, 50, 100], axis_neuron=4)
descriptor = Feature(3, neuron=[20, 40, 80], axis_neuron=4)

# Decriptors for these coordinates
#out1 = descriptor(coords, atom_types)

# Using local function
#out1 = debug_feature_forward(descriptor, coords, atom_types)


#
# debug local
#
num_batches, num_channels, _ = coords.size()
# num_batches: number of data points / coordinates
# num_channels: number of atoms

loc_env_r, loc_env_a = local_environment(coords)


# Variables in descriptor
n_types = descriptor.n_types
local_embedding = descriptor.local_embedding
axis_neuron = descriptor.axis_neuron

# make multiples?
neighbor_types = atom_types.repeat(num_channels, 1)
mask = ~torch.eye(num_channels, dtype=torch.bool, device=coords.device)
# mask is inverse identity matrix?
neighbor_types = torch.masked_select(neighbor_types, mask).view(num_channels, -1)
# Needed to satisfy permutational symmetry????

indices = ((atom_types * n_types).unsqueeze(-1) + neighbor_types).view(-1)
# unsqueeze to add singleton dimension

output, _ = local_embedding((loc_env_r.view(num_batches, -1, 1), indices))
# channel/indices is ignored in the output. It is returned again (?).

# reshape?
output = output.view(num_batches, num_channels, num_channels - 1, -1)

output = torch.transpose(output, 2, 3) @ (loc_env_a @ (torch.transpose(loc_env_a, 2, 3) @ output[..., :axis_neuron]))
output = output.view(num_batches, num_channels, -1)

out1 = output # assign

print(f"shape of out1: {out1.shape}")
print(f"sum of out1: {out1.sum()}")

