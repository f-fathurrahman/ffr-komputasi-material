import os
from utilities_database import Database

def debug_run_train(ff, TrainData=None):

    assert TrainData is not None, "TrainData can't be None for train mode."

    # Instantiate model
    print("Initializing model ...")
    ff._MODEL(ff._model)
    print("... Done initializing model")
    # ff._model is a dict

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
    #
    trainDB.close()
    
    print("=========================== Training =============================\n")
    #ff.model.train('Train_db', optimizer=ff.optimizer)
    model_nn_train(ff.model, 'Train_db', optimizer=ff.optimizer)

    #ff.model.save_checkpoint(des_info=ff._descriptors)


    #print(f"==================== Evaluating Training Set ====================\n")
    #train_stat = ff.model.evaluate('Train_db', figname='Train.png')
    train_stat = None
    test_stat = None

    return (train_stat, test_stat)


import shelve

def model_nn_train(model_nn, TrainData, optimizer):
    # NOTE: TrainData is now a name of database: Train_db

    print("Enter model_nn_train, TrainData = ", TrainData)

    """ Training of Neural Network Potential. """
    # If batch_size is None and optimizer is Adam or SGD, 
    # then batch_size equals total structures.
    if optimizer['method'] in ['sgd', 'SGD', 'Adam', 'adam', 'ADAM']:
        if model_nn.batch_size == None:
            db = shelve.open(model_nn.path+TrainData)
            model_nn.batch_size = len(db.keys())
            db.close()

    preprocess(model_nn, TrainData)
    
    # Calculate total number of parameters.
    model_nn.total_parameters = 0
    for element in model_nn.elements:
        for i, hl in enumerate(model_nn.hiddenlayers[element]):
            if i == 0:
                model_nn.total_parameters += (model_nn.no_of_descriptors+1)*hl
            else:
                model_nn.total_parameters += (model_nn.hiddenlayers[element][i-1]+1)*hl
    
    if model_nn.restart is None:
        # Creating Neural Network architectures.
        model_nn.models = {}
        for element in model_nn.elements: # Number of models depend on species
            m = 'nn.Sequential('
            for i, act in enumerate(model_nn.activation[element]):
                if i == 0:
                    m += f'nn.Linear({model_nn.no_of_descriptors}, \
                            {model_nn.hiddenlayers[element][i]}), '
                else:
                    m += f'nn.Linear({model_nn.hiddenlayers[element][i-1]}, \
                            {model_nn.hiddenlayers[element][i]}), '
                                
                if act == 'Linear':
                    continue
                else:
                    m += f'nn.{act}(), '
            m += f')'

            model_nn.models[element] = eval(m).double().to(model_nn.device)

        model_nn.regressor = Regressor(optimizer['method'], optimizer['parameters'])
        model_nn.optimizer = model_nn.regressor.regress(models=model_nn.models)

    else:
        # Look for previously saved models and continue optimizing from the last checkpoint.
        model_nn.load_checkpoint(filename=model_nn.restart, 
                                method=optimizer['method'], args=optimizer['parameters'])
            
    print(f"No of structures   : {model_nn.no_of_structures}")
    print(f"No of descriptors  : {model_nn.no_of_descriptors}")
    print(f"No of parameters   : {model_nn.total_parameters}")
    print(f"No of epochs       : {model_nn.epoch}")
    print(f"Optimizer          : {optimizer['method']}")
    print(f"Force_coefficient  : {model_nn.force_coefficient}")
    print(f"Stress_coefficient : {model_nn.stress_coefficient}")
    print(f"Batch_size         : {model_nn.batch_size}\n")


import numpy as np
import torch
from torch.utils import data
from models_neuralnetwork import Dataset
import gc

def preprocess(model_nn, TrainData):
    """ Preprocess TrainData to a convenient format for Neural Network training. """

    print("--- Enter preprocess ---")

    #if os.path.exists(model_nn.path+"drange.npy"):
    #    model_nn.drange = np.load(model_nn.path+"drange.npy", allow_pickle=True)[0]
    #else:
    #    model_nn.drange = get_descriptors_range(model_nn, TrainData)
    
    #model_nn.plot_hist(descriptors, figname=model_nn.path+"histogram.png", figsize=(12, 24))

    # force recalculate drange
    model_nn.drange = get_descriptors_range(model_nn, TrainData)

    normalize(model_nn, TrainData, model_nn.drange, model_nn.unit)

    model_nn.get_stress_group(TrainData)

    model_nn.softmax = model_nn._SOFTMAX(TrainData, beta=model_nn.softmax_beta)

    model_nn.data = data.DataLoader(Dataset(model_nn.path+TrainData, model_nn.softmax,
                                        model_nn.device, model_nn.memory, 
                                        model_nn.force_coefficient, 
                                        model_nn.stress_coefficient),
                                batch_size=model_nn.batch_size,
                                shuffle=model_nn.shuffle,
                                collate_fn=model_nn.collate_fn,)

    gc.collect()



def get_descriptors_range(model_nn, data):

    print("--- Enter get_descriptors_range")

    # force recompute
    _DRANGE = {}
    db = shelve.open(model_nn.path+data)
    no_of_structures = len(list(db.keys()))
    # keys are string

    print("data = ", data)
    print("file db = ", model_nn.path + data)
    print("no_of_structures = ", no_of_structures)

    for i in range(no_of_structures):
        dic = db[str(i)]
        for j, descriptor in enumerate(dic['x']):
            element = dic['elements'][j]
            if element not in _DRANGE.keys():
                _DRANGE[element] = np.asarray([np.asarray([__, __]) \
                                    for __ in descriptor])
            else:
                assert len(_DRANGE[element]) == len(descriptor)
                for j, des in enumerate(descriptor):
                    if des < _DRANGE[element][j][0]:
                        _DRANGE[element][j][0] = des
                    elif des > _DRANGE[element][j][1]:
                        _DRANGE[element][j][1] = des
    db.close()
    np.save(model_nn.path+'drange.npy', [_DRANGE])

    return _DRANGE

from copy import deepcopy

def normalize(model_nn, data, drange, unit, norm=[0., 1.]):

    db1 = shelve.open(model_nn.path+data)
    model_nn.no_of_structures = len(list(db1.keys()))
    model_nn.no_of_descriptors = db1['0']['x'].shape[1]
    # number of descriptors always the same for all structures, regardless of number of atoms? 

    db2 = shelve.open(model_nn.path+data+'_norm') # open file for writing

    print("no_of_structures = ", model_nn.no_of_structures)
    print("no_of_descriptors = ", model_nn.no_of_descriptors)

    for i in range(model_nn.no_of_structures):
        d = {'x': {}, 'dxdr': {}, 'seq': {}, 'rdxdr': {}}
        descriptor = db1[str(i)]

        for element in model_nn.elements:

            print("element = ", element)

            _drange = drange[element]
            scale = (norm[1] - norm[0]) / (_drange[:, 1] - _drange[:, 0])
            print("scale.shape = ", scale.shape) # this is a matrix?
            print("_drange = ", _drange.shape)
            
            res_e = np.where(np.array(descriptor['elements'])==element)
            print("res_e = ", res_e)
            e = np.where(np.array(descriptor['elements'])==element)[0]
            print("e = ", e.shape)
            try:
                ee0 = np.where(descriptor['seq'][:,0]==e[0])[0][0]
                ee1 = np.where(descriptor['seq'][:,0]==e[-1])[0][-1]
            except:
                print("Exception ...")
                ee0, ee1 = 0, 1
            
            #print("descriptor seq = ", descriptor["seq"])
            print("ee0 = ", ee0)
            print("ee1 = ", ee1)
            i_size = list(descriptor['elements']).count(element) # number of atoms in type i
            m_size = ee1 - ee0 + 1   # number of pairs
            print("m_size = ", m_size)
            print("i_size = ", i_size)

            d['x'][element] = torch.zeros([i_size, model_nn.no_of_descriptors], dtype=torch.float64)
            #d['dxdr'][element] = torch.zeros([m_size, model_nn.no_of_descriptors, 3], dtype=torch.float64) ####
            d['dxdr'][element] = np.zeros([m_size, model_nn.no_of_descriptors, 3]) ####
            d['rdxdr'][element] = torch.zeros([i_size, model_nn.no_of_descriptors, 6], dtype=torch.float64)
            d['seq'][element] = None

            if e.size > 0:
                des = norm[0] + np.einsum('j,ij->ij', scale, (descriptor['x'][e[0]:e[-1]+1] - np.expand_dims(_drange[:, 0], 0)))
                dess = np.einsum('j,ijk->ijk', scale, descriptor['rdxdr'][e[0]:e[-1]+1])
                desp = np.einsum('j,ijk->ijk', scale, descriptor['dxdr'][ee0:ee1+1])
                
                print("des.shape = ", des.shape)

                d['x'][element] += torch.from_numpy(des)
                d['seq'][element] = deepcopy(descriptor['seq'][ee0:ee1+1])
                #d['seq'][element] = torch.from_numpy(descriptor['seq'][ee0:ee1+1])

                if unit == 'eV':
                    #d['dxdr'][element] += torch.from_numpy(desp)
                    d['dxdr'][element] += desp
                    d['rdxdr'][element] += torch.from_numpy(dess)
                else:
                    d['dxdr'][element] += torch.from_numpy(0.529177 * desp)
                    d['rdxdr'][element] += torch.from_numpy(0.529177 * dess)
            #d['seq'][element][:, 0] -= torch.min(d['seq'][element][:, 0]) #adjust the initial position
                d['seq'][element][:, 0] -= np.min(d['seq'][element][:, 0]) #adjust the initial position
        
        x, dxdr, seq = d['x'], d['dxdr'], d['seq']
        n_atoms = sum(len(value) for value in x.values())
        new_dxdr = {}
        for element in x.keys():
            if x[element].nelement() > 0:
                tmp = np.zeros([len(x[element]), n_atoms, x[element].shape[1], 3])
                for _m in range(n_atoms):
                    rows = np.where(seq[element][:,1]==_m)[0]
                    tmp[seq[element][rows, 0], _m, :, :] += dxdr[element][rows, :, :]
                new_dxdr[element] = torch.from_numpy(tmp)
            else:
                new_dxdr[element] = torch.empty((30,3), dtype=torch.float)
        d['dxdr'] = new_dxdr

        db2[str(i)] = d
        
    db1.close()
    db2.close()