from ase.build import fcc100, add_adsorbate
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
#from ase.optimize import BFGS

from MyBFGS import MyBFGS
from MyMDMin import MyMDMin
from MyFIRE import MyFIRE

slab = fcc100("Al", size=(2, 2, 3))
add_adsorbate(slab, "Au", 1.7, "hollow")
slab.center(axis=2, vacuum=4.0)
slab.set_calculator(EMT())

# Fix second and third layers:
mask = [atom.tag > 1 for atom in slab]
slab.set_constraint(FixAtoms(mask=mask))

# 1.2. Optimize initial and final end-points.

# Initial end-point:
geoopt = MyBFGS(slab, trajectory="initial.traj")
#geoopt = MyMDMin(slab, trajectory="initial.traj")
#geoopt = MyFIRE(slab, trajectory="initial.traj")
geoopt.run(fmax=0.01)

print("Optimized positions: ")
print(slab.get_positions())


