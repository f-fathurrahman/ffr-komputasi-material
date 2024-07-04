import numpy as np
import pyscf.pbc.gto as pbcgto
import pyscf.pbc.scf as pbcscf
from pyscf.pbc.tools import pyscf_ase

import ase
import ase.lattice

ase_atom = ase.lattice.cubic.BodyCenteredCubic("Fe", latticeconstant=2.856)
print("Cell volume = ", ase_atom.get_volume())

# Initialize Cell object
cell = pbcgto.Cell()
cell.verbose = 5
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atoms)
cell.a = ase_atoms.cell  # we are using angstrom here
cell.basis = "gth-dzv"
cell.pseudo = "gth-pbe"
cell.nspin = 2
# need to do this because we modify cell object after call to its constructor
cell.build()

nks = [3,3,3]
mf = pbcscf.KUKS(cell, cell.make_kpts(nks))
mf.verbose = 5 # ???
mf = pbcscf.addons.smearing_(mf, sigma=0.1, method="fermi")
mf.xc = "pbe"

res1 = mf.kernel()
