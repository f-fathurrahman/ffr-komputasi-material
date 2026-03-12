from ase.io import read
import numpy as np

db = read("GABUNG01.xyz", ":")
Ndata = len(db)
Ntrain = int(0.7*Ndata)
idx_all = np.arange(Ndata)
np.random.shuffle(idx_all)
idx_train = idx_all[:Ntrain]
idx_test = idx_all[Ntrain:]

f = open("train_data.xyz", "w")
for i in idx_train:
    db[i].write("train_data.xyz", append=True)
f.close()

f = open("test_data.xyz", "w")
for i in idx_test:
    db[i].write("test_data.xyz", append=True)
f.close()