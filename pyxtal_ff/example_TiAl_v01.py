from pyxtal_ff_01 import PyXtal_FF
import os

train_data = "DATASET_others/TiAl_2_v01.xyz"

descriptor = {
    'Rc': 4.0, # make Rc smaller to speed up the calculation
    'type': 'SO3',
    'parameters': {'nmax': 4, 'lmax': 3},
    'ncpu': 1,
}
 
model = {
    'system' : ['Ti', 'Al'],
    'hiddenlayers': [30, 30],
    'force_coefficient': 1.0,
    'epoch': 100,
    'batch_size': 32,
    'path': 'TEMP_TiAl_Rc_v01/',
    'memory': 'out',
    'optimizer': {'method': 'LBFGS'}
}

#------------------------- Run NN calculation ------------------------------
ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode='train', TrainData=train_data)
