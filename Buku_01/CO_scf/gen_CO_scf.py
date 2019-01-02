from ase import Atoms
from ase.units import Bohr
import numpy as np
from qeManager.pwscf import *
from ase.build import molecule

atoms = molecule("CO")
atoms.set_pbc([True,True,True])
cell = np.array([16.0,16.0,16.0])*Bohr
atoms.set_cell(cell)
atoms.center()

# Create CONTROL namelist of PWSCF input
ctrl_NL = ControlNameList() # using default parameters
ctrl_NL.pseudo_dir = "/home/efefer/pseudo" # modify pseudo_dir directly
ctrl_NL.write_all() # write all parameters, including defaults to stdout

# also do the same for SYSTEM and ELECTRONS

sys_NL = SystemNameList(atoms)

#sys_NL.occupations = "smearing"
#sys_NL.smearing = "mv"
#sys_NL.degauss = 0.001

sys_NL.ecutwfc = 30.0
sys_NL.ecutrho = 4*sys_NL.ecutwfc
sys_NL.write_all()

elec_NL = ElectronsNameList()
elec_NL.mixing_beta = 0.1
elec_NL.write_all()

# pseudopotentials list
pspFiles = ["C_ONCV_PBE-1.0.upf", "O_ONCV_PBE-1.0.upf"]

# write ATOMIC_SPECIES card
write_atomic_species(atoms, pspFiles=pspFiles)

# write ATOMIC_POSITIONS
write_atomic_positions(atoms)

# write CELL_PARAMETERS
write_cell(atoms)

