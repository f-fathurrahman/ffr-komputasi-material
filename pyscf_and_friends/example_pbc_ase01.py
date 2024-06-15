import numpy as np
import pyscf.pbc.gto as pbcgto
import pyscf.pbc.dft as pbcdft
from pyscf.pbc.tools import pyscf_ase

import ase
import ase.lattice
from ase.lattice.cubic import Diamond

ase_atoms = Diamond(symbol="C", latticeconstant=3.5668)
print("Cell volume = ", ase_atoms.get_volume())

# Initialize Cell object
cell = pbcgto.Cell()
cell.verbose = 5
cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atoms)
cell.a = ase_atom.cell  # we are using angstrom here
cell.basis = "gth-szv"
cell.pseudo = "gth-pade"
# need to do this because we modify cell object after call to its constructor
cell.build()

mf = pbcdft.RKS(cell)
mf.xc = "lda,vwn"
# mf.xc = "pbe"  # gga,pbe did not work?

res1 = mf.kernel()
