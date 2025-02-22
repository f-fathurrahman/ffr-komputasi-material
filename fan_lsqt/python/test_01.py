import numpy as np
import matplotlib.pyplot as plt

from mod_lsqt_01 import *

# (1) Prepare some parameters
Nx = 50000 # 0                       # length
Ny = 2                            # width
W = 0                             # no disorder in this example
E_max = 3.1                       # energy scaling factor
M = 1000                          # number of Chebyshev moments
E = np.linspace(-3, 3, num = 601) # energy points
dt = np.ones(10)                  # time steps

# (2) Calculate physical quantities
DOS, VAC, sigma_from_VAC, MSD, sigma_from_MSD = lsqt(Nx, Ny, W, M, E_max, E, dt)
v_F = np.sqrt(VAC[0, :])                        # Fermi velocity
g = Ny * DOS * v_F * 0.5                        # conductance
g *= 2.0 * np.pi                                # from e^2/hbar to e^2/h


"""
# (3) Get similar plots as in Figs. 5 and 6 in the review paper

plt.figure
plt.plot(np.arange(10.0), VAC[:, 301], 'bo')
plt.xlabel('Time ($\hbar/\gamma$)')
plt.ylabel('VAC ($a^2\gamma^2/\hbar^2$)')
plt.ylim(0, 4)
plt.show()

plt.figure
plt.plot(np.arange(1.0, 11.0, 1.0), MSD[:, 301], 'bo')
plt.xlabel('Time ($\hbar/\gamma$)')
plt.ylabel('MSD ($a^2$)')
plt.show()

plt.figure
plt.plot(E, DOS)
plt.xlabel('Energy ($\gamma$)')
plt.ylabel('DOS ($1/\gamma/a^2$)')
plt.show()

plt.figure
plt.plot(E, v_F)
plt.xlabel('Energy ($\gamma$)')
plt.ylabel('v ($a\gamma/\hbar$)')
plt.show()

plt.figure
plt.plot(E, g)
plt.xlabel('Energy ($\gamma$)')
plt.ylabel('g ($e^2/h$)')
plt.show()
"""
