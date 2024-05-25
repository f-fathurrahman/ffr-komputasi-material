from my_sgdml.intf.ase_calc import SGDMLCalculator

from ase.io import read
from ase.optimize import QuasiNewton
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution, Stationary, ZeroRotation)
from ase.md.verlet import VelocityVerlet
from ase import units

model_path = "m_ethanol.npz"
calc = SGDMLCalculator(model_path)

mol = read("ethanol.xyz")
mol.set_calculator(calc)

# do a quick geometry relaxation
qn = QuasiNewton(mol)
qn.run(1e-4, 100)

# set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(mol, 300 * units.kB)
Stationary(mol) # zero linear momentum
ZeroRotation(mol) # zero angular momentum

# run MD with constant energy using the velocity verlet algorithm
dyn = VelocityVerlet(mol, 0.2 * units.fs, trajectory="md.traj")  # 0.2 fs time step.

def printenergy(a):
    # function to print the potential, kinetic and total energy
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  "
            "Etot = %.3feV" % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

# now run the dynamics
printenergy(mol)
for i in range(10000):
    dyn.run(10)
    printenergy(mol)
