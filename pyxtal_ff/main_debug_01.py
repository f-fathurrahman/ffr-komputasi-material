from my_pyxtal_ff.pyxtal_ff_01 import PyXtal_FF
import shelve

def get_desc_list(model_nn):
    db = shelve.open(model_nn.path + "Train_db")
    no_of_structures = len(list(db.keys()))
    DESC_LIST = []
    for i in range(no_of_structures):
        DESC_LIST.append(db[str(i)])
    db.close()
    return DESC_LIST


# This is extended XYZ file
#train_data = "DATASET_N2H4_v2/TEMP_ATOMS_TRAIN.xyz"
train_data = "DATASET_N2H4_v2/TEMP_2DATA.xyz"
path_model = "LOGDIR_N2H4_v2_debug_01/" # need trailing /

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

#from debug_run_train_01 import debug_run_train
#debug_run_train(ff, TrainData=train_data)


print("Initializing model ...")
ff._MODEL(ff._model)
print("... Done initializing model")
# ff._model is a dict

TrainData = train_data # set argument

import os
from my_pyxtal_ff.utilities_database import Database

# Calculate descriptors.
ff._descriptors.update({'N': ff._descriptors['N_train']})
if not os.path.exists(ff.path+'Train_db.dat') and not os.path.exists(ff.path+'Train_db.db'):
    trainDB = Database(name=ff.path+'Train_db')
    trainDB.store(TrainData, ff._descriptors, True, ff.path+'ase.db')
else:
    # The database is already exist, just load it
    # XXX In several case this will error if the database is not complete
    trainDB = Database(name=ff.path+'Train_db')
    trainDB.store(TrainData, ff._descriptors, False)
trainDB.close()

#from debug_run_train_01 import model_nn_train
#model_nn_train(ff.model, 'Train_db', optimizer=ff.optimizer)

from debug_run_train_01 import get_descriptors_range, normalize
ff.model.drange = get_descriptors_range(ff.model, "Train_db")
normalize(ff.model, "Train_db", ff.model.drange, ff.model.unit)

