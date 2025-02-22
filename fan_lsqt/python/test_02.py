import numpy as np
import matplotlib.pyplot as plt

from mod_lsqt_01 import *


# (1) Prepare some parameters
Nx = 20000                         # length
Ny = 50                            # width
W = 1.0                            # disordered
E_max = 5.1                        # energy scaling factor
M = 1000                           # number of Chebyshev moments
E = np.linspace(-5, 5, num = 1001) # energy points
dt = np.ones(20) * 2.0             # time steps

# (2) Calculate physical quantities
DOS, VAC, sigma_from_VAC, MSD, sigma_from_MSD = lsqt(Nx, Ny, W, M, E_max, E, dt)

# (3) Get plots similar to those in Fig. 7 in the review paper
plt.figure
plt.plot(np.arange(0.0, 40.0, 2.0), VAC[:, 501], 'bo')
plt.xlabel('Time ($\hbar/\gamma$)')
plt.ylabel('VAC ($a^2\gamma^2/\hbar^2$)')
plt.show()

plt.figure
plt.plot(np.arange(2.0, 42.0, 2.0), MSD[:, 501], 'rd')
plt.xlabel('Time ($\hbar/\gamma$)')
plt.ylabel('MSD ($a^2$)')
plt.show()

plt.figure
plt.plot(np.arange(0.0, 40.0, 2.0), sigma_from_VAC[:, 501], 'bo')
plt.plot(np.arange(1.0, 41.0, 2.0), sigma_from_MSD[:, 501], 'rd')


plt.xlabel('Time ($\hbar/\gamma$)')
plt.ylabel('$\sigma$ ($e^2/h$)')
plt.show()

plt.figure
plt.plot(E, np.amax(sigma_from_VAC, axis=0), 'bo')
plt.plot(E, np.amax(sigma_from_MSD, axis=0), 'rd')
plt.xlabel('Energy ($\gamma$)')
plt.ylabel('$\sigma_{sc}$ ($e^2/h$)')
plt.show()




test_lsqt_v01()