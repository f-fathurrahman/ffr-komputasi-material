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
    ff.model.save_checkpoint(des_info=ff._descriptors)
        
    print(f"==================== Evaluating Training Set ====================\n")
    train_stat = ff.model.evaluate('Train_db', figname='Train.png')
    
    test_stat = None

    return (train_stat, test_stat)