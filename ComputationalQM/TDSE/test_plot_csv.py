import numpy as np
import sys

filename = sys.argv[1]
prefix = filename.replace(".csv", "")

dat = np.loadtxt(filename, delimiter=";")

Nsamples = dat.shape[1] - 1

import matplotlib.pyplot as plt

t = dat[:,0]
plt.clf()
for i in range(1,Nsamples+1):
    plt.plot(t, dat[:,i])
plt.savefig(prefix + ".png", dpi=150)
