import sys
import numpy as np

import numpy.random
numpy.random.seed(1234)

from my_sgdml.train import GDMLTrain

dataset = np.load("DATASET/ethanol_dft.npz")
n_train = 200

gdml_train = GDMLTrain(max_processes=1)
task = gdml_train.create_task(dataset, n_train,\
        valid_dataset=dataset, n_valid=1000,\
        sig=20, lam=1e-10, use_sym=False)

model = gdml_train.train(task)

#try:
#except Exception, err:
#    sys.exit(err)
#else:
#    np.savez_compressed('m_ethanol.npz', **model)

