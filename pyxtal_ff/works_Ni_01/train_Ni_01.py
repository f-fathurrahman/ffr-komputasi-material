import sys
sys.path.append('../')

from my_pyxtal_ff.pyxtal_ff_01 import PyXtal_FF

train_data = "ALL_ATOMS_T500.xyz"
path_model = "LOGDIR_Ni_fcc/" # need trailing /

descriptor = {
    "Rc": 4.0, # make Rc smaller to speed up the calculation
    "type": "SO3",
    "parameters": {"nmax": 4, "lmax": 3},
    "ncpu": 1,
}

model = {
    "algorithm": "PR",
    "system" : ["Ni"],
    "path": path_model,
    'force_coefficient': 0.001,
}

""""
model = {
    "system" : ["Ni"],
    "hiddenlayers": [10],
    "force_coefficient": 0.00001,
    "epoch": 50,
    "batch_size": 32,
    "path": path_model,
    "memory": "out",
    "optimizer": {"method": "LBFGS"}
}
"""

ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode="train", TrainData=train_data)
