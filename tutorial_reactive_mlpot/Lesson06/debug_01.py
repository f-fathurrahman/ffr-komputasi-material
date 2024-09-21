import numpy as np
import torch
from my_deeppot import *
from prepare_data import *

## XXX atom_types_all is not used in any of Feature, Fitting, or DeepPot instances

qm_coord, atom_types_all, energy, qm_gradient = prepare_data()

descriptor = Feature(3, neuron=[25, 50, 100], axis_neuron=4)
fitting_net = Fitting(3, descriptor.output_length)
model = DeepPot(descriptor, fitting_net, learning_rate=5e-4)

# Evaluate for one data point
#with torch.no_grad():
#ene_pred, grad_pred = model(qm_coord[0:1,:,:], atom_types_all[0])
#print(f"ene_pred = {ene_pred}")

#
# Prepare arguments
#

# coordinate, only one data point
coords = qm_coord[0:1,:,:]
atom_types = atom_types_all[0]
coords.requires_grad_()

# Decriptors for these coordinates
descriptors = descriptor(coords, atom_types)

# forward for fitting_net
atomic_energies = fitting_net( (descriptors, atom_types) )
print(f"atomic_energies.shape = {atomic_energies.shape}")

energy = torch.unbind(torch.sum(atomic_energies, dim=1))

gradient, = torch.autograd.grad(energy, [coords], create_graph=True)
ene_pred = torch.hstack(energy)
grad_pred = gradient

print(f"ene_pred = {ene_pred}")