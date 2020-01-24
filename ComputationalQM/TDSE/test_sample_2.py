import numpy as np
import matplotlib.pyplot as plt

from math import exp

import matplotlib.style
matplotlib.style.use("seaborn")

from scipy.stats import rv_continuous

def psi_0(x):
    vx = 10     # value of the initial velocity
    f = 0.j + np.exp(-((x+15)/4)**2)  # Gaussian profile
    f = f*np.exp(1j*vx*x)   # Multiply by an x-dependent phase to introduce velocity
    return f

class PDF_psi0_gen(rv_continuous):

    def _pdf(self, x):
        #return np.exp(-(x+15)**2 / 2.) / np.sqrt(2.0 * np.pi)
        f = np.exp( -(x + 15.0)**2 /2.0 )/np.sqrt(2.0 * np.pi)
        return f #**2

initial_gen = PDF_psi0_gen(name="initial")
x = initial_gen.rvs(size=100)
plt.clf()
plt.plot(x, np.zeros(100), marker="o")
plt.hist(x, bins=30)
plt.savefig("IMG_initial_distrib.png", dpi=150)

plt.clf()
x = np.linspace(-30.0, -5, 1000)
y = abs(psi_0(x))**2
plt.plot(x, y)
plt.savefig("IMG_psi0_2.png", dpi=150)

