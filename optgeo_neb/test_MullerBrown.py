import copy

from ase import Atoms
from MullerBrown import MullerBrown
from ase.optimize import BFGS, MDMin

ase_calculator = MullerBrown()

# "Molecular" structures:
# Using one atom only, position at xy plane only (with z=0)
initial_structure = Atoms('C', positions=[(-0.55, 1.30, 0.0)])
final_structure = Atoms('C', positions=[(0.626, 0.025, 0.0)])

# Why use deepcopy here?
initial_structure.set_calculator(copy.deepcopy(ase_calculator))
final_structure.set_calculator(copy.deepcopy(ase_calculator))

# Test calculating potential energy and forces
print("Potential energy = ", initial_structure.get_potential_energy())
print("Forces           = ", initial_structure.get_forces())


print("Optimization\n")
print("initial_structure positions: ", initial_structure.get_positions())
initial_opt = BFGS(initial_structure, trajectory='initial_optimized.traj')
initial_opt.run(fmax=0.01)
print("initial_structure positions: ", initial_structure.get_positions())
#print("Pass here")

print()
print("final_structure positions: ", final_structure.get_positions())
final_opt = BFGS(final_structure, trajectory='final_optimized.traj')
final_opt.run(fmax=0.01)
print("final_structure positions: ", final_structure.get_positions())
