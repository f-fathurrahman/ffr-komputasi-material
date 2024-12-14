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
    ff.model.train('Train_db', optimizer=ff.optimizer)
    model_nn_train(ff.model, 'Train_db', optimizer=ff.optimizer)

    #ff.model.save_checkpoint(des_info=ff._descriptors)


    #print(f"==================== Evaluating Training Set ====================\n")
    #train_stat = ff.model.evaluate('Train_db', figname='Train.png')
    train_stat = None
    test_stat = None

    return (train_stat, test_stat)


import shelve

def model_nn_train(model_nn, TrainData, optimizer):

    print("Enter model_nn_train")

    """ Training of Neural Network Potential. """
    # If batch_size is None and optimizer is Adam or SGD, 
    # then batch_size equals total structures.
    if optimizer['method'] in ['sgd', 'SGD', 'Adam', 'adam', 'ADAM']:
        if model_nn.batch_size == None:
            db = shelve.open(model_nn.path+TrainData)
            model_nn.batch_size = len(db.keys())
            db.close()

    model_nn.preprocess(TrainData)
    
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

