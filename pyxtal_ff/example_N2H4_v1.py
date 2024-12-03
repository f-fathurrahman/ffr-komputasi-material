from pyxtal_ff_01 import PyXtal_FF
import os

# This is a VASP output file
train_data = "DATASET_N2H4_v1/TEMP_ATOMS_TRAIN.xyz"

descriptor = {
    'Rc': 3.0, # make Rc smaller to speed up the calculation
    'type': 'SO3',
    'parameters': {'nmax': 4, 'lmax': 3},
    'ncpu': 1,
}
 
model = {
    'system' : ['N', 'H'],
    'hiddenlayers': [30, 30],
    'force_coefficient': 1.0,
    'epoch': 100,
    'batch_size': 32,
    'path': 'TEMP_N2H4_database_v2/',
    'memory': 'out',
    'optimizer': {'method': 'LBFGS'}
}

#------------------------- Run NN calculation ------------------------------
ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode='train', TrainData=train_data)
