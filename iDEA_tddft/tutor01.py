import numpy as np
import matplotlib.pyplot as plt

import iDEA

# A predefined system
#atom = iDEA.system.systems.atom # too large, 4 electrons

#__x2 = np.linspace(-20, 200, 200)
#atom = iDEA.system.System(
#    __x2, -2.0 / (abs(__x2) + 1.0), iDEA.interactions.softened_interaction(__x2), "ud"
#)

# Harmonic potential
x = np.linspace(-10, 10, 150)
v_ext = 0.5 * 0.25**2 * x**2
v_int = iDEA.interactions.softened_interaction(x)
atom = iDEA.system.System(x, v_ext, v_int, electrons='ud') # also test ud

#ground_state = iDEA.methods.interacting.solve(atom, k=0)
#print("Total energy = ", ground_state.energy)
#print(ground_state.space.shape)
#print(ground_state.spin.shape)
#print(ground_state.full.shape)


current_method = iDEA.methods.hartree_fock
ground_state = current_method.solve(atom, k=0)
Etot = current_method.total_energy(atom, ground_state)
print("Total energy (Hartree-Fock) = ", Etot)

#current_method = iDEA.methods.hartree
#ground_state = current_method.solve(atom, k=0)
#Etot = current_method.total_energy(atom, ground_state)
#print("Total energy (Hartree) = ", Etot)

current_method = iDEA.methods.lda
ground_state = current_method.solve(atom, k=0)
Etot = current_method.total_energy(atom, ground_state)
print("Total energy (LDA) = ", Etot)

#vmax = np.max(ground_state.space.real)
#vmin = -np.max(ground_state.space.real)
#plt.clf()
#plt.imshow(ground_state.space.real, cmap="seismic")
#plt.show()