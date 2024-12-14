from pyxtal_ff_01 import PyXtal_FF

# This is extended XYZ file
train_data = "DATASET_N2H4_v2/TEMP_ATOMS_TRAIN.xyz"
path_model = "TEMP_N2H4_v2_debug_01/" # need trailing /

descriptor = {
    'Rc': 4.0, # make Rc smaller to speed up the calculation
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
    'path': path_model,
    'memory': 'out',
    'optimizer': {'method': 'LBFGS'}
}

ff = PyXtal_FF(descriptors=descriptor, model=model)
#ff.run(mode='train', TrainData=train_data)

from debug_run_train_01 import debug_run_train
debug_run_train(ff, TrainData=train_data)
