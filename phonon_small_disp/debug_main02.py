from ase import Atoms, Atom
from my_emt import MyEMT
from my_phonons import MyPhonons

# Setup crystal and EMT calculator
atoms = Atoms()
atoms.append( Atom("Al", [0.0, 0.0, 0.0]) )
atoms.set_cell([1.0, 20.0, 20.0])
atoms.set_pbc([True, True, True])

N = 8
ph = MyPhonons(atoms, MyEMT(), supercell=(N, 1, 1), delta=0.05)
ph.run()

# Read forces and assemble the dynamical matrix
from m_debug_phonon import *
debug_phonon_read(ph, acoustic=True)

print("ph.D_N.shape = ", ph.D_N.shape)

ph.clean()

# Bandpath
# special_points={AGMRXZ}
path = atoms.cell.bandpath("GZ", npoints=100)
bs = ph.get_band_structure(path)

# Plot the band structure and DOS:
import matplotlib.pyplot as plt
plt.style.use("dark_background")

fig = plt.figure(1, figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.035
bs.plot(ax=ax, emin=0.0, emax=0.005)


"""
dos = ph.get_dos(kpts=(20, 1, 1)).sample_grid(npts=100, width=1e-3)
dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color="grey",
                   edgecolor="k", lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)
"""

fig.savefig("TEMP_Al_1d_phonon.png", dpi=150)
