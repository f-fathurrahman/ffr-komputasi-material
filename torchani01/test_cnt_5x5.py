from ase.build import nanotube
from ase.md.langevin import Langevin
from ase.optimize import BFGS
from ase import units

import torch
import my_torchani

# Now let's set up a crystal
atoms = nanotube(5, 5, length=10, bond=1.420, symbol='C')
print(len(atoms), "atoms in the cell")
atoms.set_pbc([True,True,True])
Lz = atoms.cell.array[2,2]
Lx = 20.0
Ly = 20.0
atoms.set_cell([Lx, Ly, Lz])
atoms.center()
atoms.write("STRUCT_cnt_5x5.xsf")

#device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = "cpu"
print("device = ", device)

# Now let's create a calculator from builtin models:
calculator = my_torchani.models.ANI1ccx().to(device).ase()

# Now let's set the calculator for ``atoms``:
atoms.set_calculator(calculator)

# Now let's minimize the structure:
print("Begin minimizing...")
opt = BFGS(atoms)
opt.run(fmax=0.001)
print()


# Now create a callback function that print interesting physical quantities:
def printenergy(a=atoms):
    # Function to print the potential, kinetic and total energy.
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))


# We want to run MD with constant energy using the Langevin algorithm
# with a time step of 1 fs, the temperature 300K and the friction
# coefficient to 0.02 atomic units.
dyn = Langevin(atoms, 1 * units.fs, 300 * units.kB, 0.2)
dyn.attach(printenergy, interval=1)

from ase.io.trajectory import Trajectory
traj = Trajectory('cnt_5x5.traj', 'w', atoms)
dyn.attach(traj.write, interval=1)

# Now run the dynamics:
print("Beginning dynamics...")
printenergy()
dyn.run(500)

