from pyxtal_ff_01 import PyXtal_FF
import os

train_data = "DATASET_others/TiAl_gabung.xyz"
path_model = "LOGDIR_TiAl_PR_sorted/" # need trailing /

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
