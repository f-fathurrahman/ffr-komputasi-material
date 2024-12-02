import pickle
import time

import numpy as np

from descriptors_SO3 import SO3


with open("DATASET_N2H4_v1/TEMP_ATOMS_TRAIN.pkl", "rb") as f:
    atoms_train = pickle.load(f)

lmax = 2
nmax = 1
rcut = 3.0
alpha = 2.0

f = SO3(nmax=nmax, lmax=lmax, rcut=rcut, alpha=alpha, derivative=True, stress=False, cutoff_function='cosine')

start1 = time.time()
for atoms in atoms_train:
    x = f.calculate(atoms)
start2 = time.time()

#print('x', x['x'])
#print('dxdr', x['dxdr'])
print('calculation time {}'.format(start2-start1))
