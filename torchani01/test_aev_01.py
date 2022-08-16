import os
import sys
import torch
import my_torchani
import numpy as np

def load_data():
    # ani_gdb_s01.h5  ani_gdb_s02.h5  ani_gdb_s03.h5  ani_gdb_s04.h5
    #dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani1-up_to_gdb4/ani_gdb_s01.h5"
    dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani-1x/sample.h5"
    dataset = my_torchani.data.load(dspath, additional_properties=("forces",))

    energy_shifter = my_torchani.utils.EnergyShifter(None)
    species_order = ["H", "C", "N", "O"]
    dataset = dataset.subtract_self_energies(energy_shifter, species_order).species_to_indices(species_order)

    DATASET_LIST2 = []
    for i,v in enumerate(dataset):
        DATASET_LIST2.append(v)

    return DATASET_LIST2


#
# Create AEV computer
#
device = "cpu"

Rcr = 5.2000e+00
Rca = 3.5000e+00

# Radial
EtaR = torch.tensor([1.6000000e+01], device=device)
ShfR = torch.tensor([
    9.0000000e-01, 1.1687500e+00, 1.4375000e+00, 1.7062500e+00,
    1.9750000e+00, 2.2437500e+00, 2.5125000e+00, 2.7812500e+00,
    3.0500000e+00, 3.3187500e+00, 3.5875000e+00, 3.8562500e+00,
    4.1250000e+00, 4.3937500e+00, 4.6625000e+00, 4.9312500e+00], device=device)

# Angular
Zeta = torch.tensor([3.2000000e+01], device=device)
ShfZ = torch.tensor([
    1.9634954e-01, 5.8904862e-01, 9.8174770e-01, 1.3744468e+00,
    1.7671459e+00, 2.1598449e+00, 2.5525440e+00, 2.9452431e+00], device=device)

EtaA = torch.tensor([8.0000000e+00], device=device)
ShfA = torch.tensor([9.0000000e-01, 1.5500000e+00, 2.2000000e+00, 2.8500000e+00], device=device)

species_order = ["H", "C", "N", "O"]
num_species = len(species_order)

aev_computer = my_torchani.AEVComputer(Rcr, Rca, EtaR, ShfR, EtaA, Zeta, ShfA, ShfZ, num_species)
aev_dim = aev_computer.aev_length
print("aev_dim = ", aev_dim)

#
# Test AEV on dataset
#
DATASET_DICTS = load_data() # load dataset as list of dict
for i in range(5000,5010):
    # Convert to tensor and resize to match AEV expectation
    sp1 = torch.tensor(
        np.expand_dims(DATASET_DICTS[i]["species"], axis=0)
    )
    coord1 = torch.tensor(
        np.expand_dims(DATASET_DICTS[i]["coordinates"], axis=0)
    )
    res = aev_computer( (sp1, coord1) ) # call forward method
    print()
    print(res.species)
    print(res.species.shape)
    print(res.aevs.shape)

