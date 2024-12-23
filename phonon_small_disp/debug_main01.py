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
from m_debug_phonon import *
debug_phonon_read(ph, acoustic=True)
ph.clean()

print("ph.D_N.shape = ", ph.D_N.shape)

path = atoms.cell.bandpath("GXULGK", npoints=100)
bs = get_band_structure(ph, path) #XXX cannot compute modes?
