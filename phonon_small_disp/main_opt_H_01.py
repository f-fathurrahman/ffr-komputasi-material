from ase import Atoms, Atom
from my_emt import MyEMT
import numpy as np

# Setup crystal and EMT calculator
atoms = Atoms()
atoms.append(Atom("H", [0.0, 0.0, 0.0]))
atoms.set_pbc([True, True, True])
atoms.calc = MyEMT()

acell_list = np.linspace(0.5, 2.0, 20)
energies = []
for acell in acell_list:
    atoms.set_cell([acell, 20.0, 20.0])
    E = atoms.get_potential_energy()
    energies.append(E)

import matplotlib.pyplot as plt
plt.plot(acell_list, energies)
plt.grid()
plt.show()
