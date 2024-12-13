from pyxtal_ff_01 import PyXtal_FF

# This is extended XYZ file
train_data = "DATASET_N2H4_v2/TEMP_ATOMS_TRAIN.xyz"
path_model = "TEMP_N2H4_v2_gabung_v3/" # need trailing /

descriptor = {
    'Rc': 4.0, # make Rc smaller to speed up the calculation
    'type': 'SO3',
    'parameters': {'nmax': 4, 'lmax': 3},
    'ncpu': 1,
}
 
model = {
    "algorithm": "PR",
    "system" : ["N", "H"],
    "path": path_model,
    'force_coefficient': 0.5,
}

#------------------------- Run NN calculation ------------------------------
ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode='train', TrainData=train_data)
