import numpy as np
import matplotlib.pyplot as plt
from excor import ExchangeCorrelation

# Some parameters
RmaxAtom = 10.0  # The end of the radial mesh (maximum r)
Nratom = 3001    # Number of points in radial mesh

XC = ExchangeCorrelation(3)

# Radial mesh, equally spaced
R0 = np.linspace(1e-10, RmaxAtom, Nratom) 
# Inverse radial mesh
Ra = R0[::-1]

Veff = -np.ones(len(Ra), dtype=float)/Ra

plt.clf()
plt.plot(R0, Veff)
plt.savefig("IMG_atom_Veff.pdf")
