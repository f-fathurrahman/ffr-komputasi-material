from pyxtal_ff_01 import PyXtal_FF
import os

#train_data = "DATASET_others/TiAl_2atoms_v01.xyz"
train_data = "DATASET_others/TiAl_gabung.xyz"
path_model = "TEMP_TiAl_gabung_v01/" # need trailing /

descriptor = {
    "Rc": 4.0, # make Rc smaller to speed up the calculation
    "type": "SO3",
    "parameters": {"nmax": 4, "lmax": 3},
    "ncpu": 1,
}


"""
symmetry = {
    "G2": {
        "eta": [0.035709, 0.071418, 0.178545, 0.35709, 0.71418, 1.78545],
        "Rs": [0]
    },
    "G4": {
        "lambda": [-1, 1],
        "zeta": [1],
        "eta": [0.035709, 0.071418, 0.178545, 0.35709]
    }
}

descriptor = {
    "type": "BehlerParrinello",
    "parameters": symmetry,
    "Rc": 5.0,
    "ncpu": 1,
}
"""


model = {
    "system" : ["Ti", "Al"],
    "hiddenlayers": [30, 30],
    "force_coefficient": 1.0,
    "epoch": 100,
    "batch_size": 32,
    "path": path_model,
    "memory": "out",
    "optimizer": {"method": "LBFGS"}
}

ff = PyXtal_FF(descriptors=descriptor, model=model)

#ff.run(mode="train", TrainData=train_data)
from debug_run_train_01 import debug_run_train
debug_run_train(ff, TrainData=train_data)