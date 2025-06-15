from ase import Atoms, Atom
from my_emt import MyEMT
from my_phonons import MyPhonons

# Simulating periodic one dimension crystal

# Setup crystal and EMT calculator
atoms = Atoms()
atoms.append(Atom("H", [0.0, 0.0, 0.0], mass=1.0))
atoms.append(Atom("H", [0.91, 0.0, 0.0], mass=1.0))
atoms.set_pbc([True, True, True])
atoms.set_cell([1.8, 20.0, 20.0])


# Probably need to optimize the geometry or unit cell?

N = 9
ph = MyPhonons(atoms, MyEMT(), supercell=(N, 1, 1), delta=0.01)
ph.run()

# Read forces and assemble the dynamical  matrix
from m_debug_phonon import *
debug_phonon_read(ph, acoustic=False, symmetrize=0)
print("ph.D_N.shape = ", ph.D_N.shape) # (N,3*Natoms,3*Natoms)
ph.clean()

# Bandpath
# special_points={AGMRXZ}
path = atoms.cell.bandpath("GZ", npoints=100)
omega_kl = band_structure(ph, path.kpts, modes=False, born=False)
bs = ph.get_band_structure(path)

import matplotlib.pyplot as plt
qpts_x = path.kpts[:,0]
Nbands = omega_kl.shape[1]
for ibnd in range(Nbands):
    all_pos = np.all(omega_kl[:,ibnd] >= 0)
    if all_pos:
        plt.plot(qpts_x, omega_kl[:,ibnd], label="ibnd-"+str(ibnd))
plt.grid()
plt.legend()
plt.show()



"""
# Plot the band structure and DOS:
import matplotlib.pyplot as plt
plt.style.use("dark_background")

# Create figure with given size
fig = plt.figure(1, figsize=(7, 4))
# Create an axis
ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])

emax = 0.035
bs.plot(ax=ax, emin=0.0, emax=0.005)


dos = ph.get_dos(kpts=(20, 1, 1)).sample_grid(npts=100, width=1e-3)
dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color="grey",
                   edgecolor="k", lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)


fig.savefig("TEMP_Al_1d_phonon.png", dpi=150)
"""