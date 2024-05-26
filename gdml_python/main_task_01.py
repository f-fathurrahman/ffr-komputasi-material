import sys
import numpy as np
from my_sgdml.train import GDMLTrain

dataset = np.load("DATASET/ethanol_dft.npz")
n_train = 200

gdml_train = GDMLTrain()
task = gdml_train.create_task(
    dataset, n_train, use_sym=True,
    valid_dataset=dataset, n_valid=1000,
    sig=20, lam=1e-10
)


R = task["R_train"]
n_train, n_atoms, _ = R.shape

import scipy.spatial.distance

i = 0
r = np.squeeze(R[i, :, :])
adj = scipy.spatial.distance.pdist(r, 'euclidean')
sadj = scipy.spatial.distance.squareform(adj)
w, v = np.linalg.eig(sadj)

# From large to smallest
v = v[:, w.argsort()[::-1]]
