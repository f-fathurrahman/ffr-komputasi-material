from ase import Atoms
from ase.units import Bohr
import numpy as np
from qeManager import ConvergenceTest, PWSCFInput
from ase.build import molecule

atoms = molecule("LiH")
atoms.set_pbc([True,True,True])
cell = np.array([16.0,16.0,16.0])*Bohr
atoms.set_cell(cell)
atoms.center()

pspFiles = ["Li.pbe-s-kjpaw_psl.0.2.1.UPF", "H.pbe-kjpaw.UPF"]
pwinput = PWSCFInput(
    atoms, pspFiles,
    filename="PWINPUT", move_atoms=False,
    gamma_only=True )

pwinput.CONTROL.pseudo_dir = "/home/efefer/pseudo"
pwinput.CONTROL.disk_io = "none"
pwinput.ELECTRONS.mixing_beta = 0.5

conv_test = ConvergenceTest(
    pwinput, what="ecutwfc",
    values=np.arange(30.0,80.0,5.0) )

#conv_test.run()

ecutrho, energies, total_forces = conv_test.read()

Ndata = len(ecutrho)
Eref = energies[-1]  # take the most converged data (the last one) as the reference
ForceRef = total_forces[-1]

print('\nConvergence in total energy:')
for i in range(Ndata):
    print("%10.5f %18.10f %18.10f" % (ecutrho[i], energies[i], energies[i]-Eref))

print('\nConvergence in total force:')
for i in range(Ndata):
    print("%10.5f %18.10f %18.10f" % (ecutrho[i], total_forces[i],
        total_forces[i]-ForceRef))
