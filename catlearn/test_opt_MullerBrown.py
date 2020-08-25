import copy

from ase import Atoms
from MullerBrown import MullerBrown
from MyBFGS import MyBFGS
from MyMDMin import MyMDMin
from MyFIRE import MyFIRE

calc = MullerBrown()

# "Molecular" structures:
# Using one atom only, position at xy plane only (with z=0)
atoms = Atoms("C", positions=[(-0.55, 1.30, 0.0)])
atoms.set_calculator(copy.deepcopy(calc))

# Test calculating potential energy and forces
print("Potential energy = ", atoms.get_potential_energy())
print("Forces           = ", atoms.get_forces())

print("Optimization\n")
print("Initial positions: ", atoms.get_positions())
#geoopt = MyBFGS(atoms, trajectory='geoopt.traj')
#geoopt = MyMDMin(atoms, trajectory='geoopt.traj')
geoopt = MyFIRE(atoms, trajectory='geoopt.traj', force_consistent=False)
geoopt.run(fmax=0.01)
print("Optimized positions: ", atoms.get_positions())
