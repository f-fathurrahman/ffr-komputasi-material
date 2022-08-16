import os
import sys
import torch
import my_torchani

torch.manual_seed(1234) # for reproducibility

# ani_gdb_s01.h5  ani_gdb_s02.h5  ani_gdb_s03.h5  ani_gdb_s04.h5
#dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani1-up_to_gdb4/ani_gdb_s01.h5"
dspath = "/home/efefer/WORKS/TORCHANI/torchani-master/dataset/ani-1x/sample.h5"
dataset = my_torchani.data.load(dspath, additional_properties=('forces',))

print("type(dataset) = ", type(dataset)) # TransformableIterable

DATASET_LIST = []
for i,v in enumerate(dataset):
    DATASET_LIST.append(v)
print("Length DATASET_LIST = ", len(DATASET_LIST))
print("Sizeof DATASET_LIST = ", sys.getsizeof(DATASET_LIST))

# Now DATASET_LIST contains list of dict
print("Keys of dict in DATASET_LIST:")
print(DATASET_LIST[0].keys())


#

energy_shifter = my_torchani.utils.EnergyShifter(None)
species_order = ["H", "C", "N", "O"]
dataset = dataset.subtract_self_energies(energy_shifter, species_order).species_to_indices(species_order)

DATASET_LIST2 = []
for i,v in enumerate(dataset):
    DATASET_LIST2.append(v)

print("Ndata DATASET_LIST2  = ", len(DATASET_LIST2))
print("Sizeof DATASET_LIST2 = ", sys.getsizeof(DATASET_LIST2))

import numpy as np
# Comparison
idx = 1000
sae = energy_shifter.self_energies
ene_orig = DATASET_LIST[idx]["energies"]
ene_shifted = DATASET_LIST2[idx]["energies"]
print("Original: ", ene_orig)
print("Shifted: ", ene_shifted)
s0 = np.sum(DATASET_LIST2[idx]["species"] == 0)
s1 = np.sum(DATASET_LIST2[idx]["species"] == 1)
s2 = np.sum(DATASET_LIST2[idx]["species"] == 2)
s3 = np.sum(DATASET_LIST2[idx]["species"] == 3)
print("ene_shifted2 = ", ene_orig - s0*sae[0] - s1*sae[1] - s2*sae[2] - s3*sae[3])

#
#training, validation = my_torchani.data.load(dspath).subtract_self_energies(energy_shifter, species_order).species_to_indices(species_order).shuffle().split(0.8, None)

#DATASET_TRAINING = []
#for i,v in enumerate(training):
#    DATASET_TRAINING.append(v)

#DATASET_VAL = []
#for i,v in enumerate(validation):
#    DATASET_VAL.append(v)
