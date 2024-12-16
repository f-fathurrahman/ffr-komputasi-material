from pyxtal_ff_01 import PyXtal_FF
import os

# This is extended XYZ file
train_data = "DATASET_N2H4_v2/TEMP_ATOMS_TRAIN.xyz"
path_model = "LOGDIR_N2H4_gabung_v1/" # need trailing /

descriptor = {
    'Rc': 4.0, # make Rc smaller to speed up the calculation
    'type': 'SO3',
    'parameters': {'nmax': 4, 'lmax': 3},
    'ncpu': 1,
}
 
model = {
    "system" : ["N", "H"],
    "hiddenlayers": [30, 30],
    "force_coefficient": 0.001,
    "epoch": 50,
    "batch_size": 32,
    "path": path_model,
    "memory": "out",
    "optimizer": {"method": "LBFGS"}
}

#------------------------- Run NN calculation ------------------------------
ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode='train', TrainData=train_data)
