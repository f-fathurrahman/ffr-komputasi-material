from ase.build import bulk

#from ase.calculators.emt import EMT
from my_emt import MyEMT

#from ase.phonons import Phonons
from my_phonons import MyPhonons

# Setup crystal and EMT calculator
atoms = bulk("Al", "fcc", a=4.05)

# Phonon calculator
N = 7
ph = MyPhonons(atoms, MyEMT(), supercell=(N, N, N), delta=0.05)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
ph.clean()

path = atoms.cell.bandpath("GXULGK", npoints=100)
bs = ph.get_band_structure(path)

dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

# Plot the band structure and DOS:
import matplotlib.pyplot as plt
plt.style.use("dark_background")

fig = plt.figure(1, figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.035
bs.plot(ax=ax, emin=0.0, emax=emax)

dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color="grey",
                   edgecolor="k", lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)

fig.savefig("TEMP_Al_phonon.png", dpi=150)
