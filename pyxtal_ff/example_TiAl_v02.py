from pyxtal_ff_01 import PyXtal_FF
import os

#train_data = "DATASET_others/TiAl_2atoms_v01.xyz"
train_data = "DATASET_others/TiAl_gabung.xyz"
path_model = "TEMP_TiAl_gabung_v02/" # need trailing /

descriptor = {
    "Rc": 4.0, # make Rc smaller to speed up the calculation
    "type": "SO3",
    "parameters": {"nmax": 4, "lmax": 3},
    "ncpu": 1,
}

model = {
    "algorithm": "PR",
    "system" : ["Ti", "Al"],
    "path": path_model,
}

ff = PyXtal_FF(descriptors=descriptor, model=model)

ff.run(mode="train", TrainData=train_data)
