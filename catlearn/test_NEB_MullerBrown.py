import copy

import ase.io
from MullerBrown import MullerBrown
from ase.optimize import BFGS, MDMin
#from ase.neb import NEB
from MyNEB import MyNEB
from MyMDMin import MyMDMin
from MyBFGS import MyBFGS

ase_calculator = MullerBrown()

n_images = 12

# Run test_MullerBrown.py to produce these files
initial_ase = ase.io.read("initial_optimized.traj")
final_ase = ase.io.read("final_optimized.traj")

images_ase = [initial_ase]

for i in range(1, n_images-1):
    image_ase = initial_ase.copy()
    image_ase.set_calculator(copy.deepcopy(ase_calculator))
    images_ase.append(image_ase)
images_ase.append(final_ase)

neb_ase = MyNEB(images_ase, climb=True, method="aseneb") # other methods are not working
#neb_ase.interpolate(method="idpp")
neb_ase.interpolate(method="linear") # default

print("Initial potential energy = ", neb_ase.get_potential_energy())
print("Initial forces")
print(neb_ase.get_forces())

#print(dir(neb_ase))
for i,image in enumerate(neb_ase.images):
    print("%3d %18.10f" % (i, image.get_potential_energy()))

#qn_ase = MDMin(neb_ase, trajectory="neb_ase.traj")
#qn_ase.run(fmax=0.05)
#print("\nSummary of the results: \n")
#atoms_ase = ase.io.read("neb_ase.traj", ":")
#n_eval_ase = int(len(atoms_ase) - 2 * (len(atoms_ase)/n_images))
#print("Number of function evaluations CI-NEB implemented in ASE:", n_eval_ase)

#bfgs_ase = BFGS(neb_ase, trajectory="neb_ase_bfgs.traj")

bfgs_ase = MyBFGS(neb_ase, trajectory="neb_ase_bfgs.traj")
bfgs_ase.run(fmax=0.05)
print("\nSummary of the results: \n")
atoms_ase = ase.io.read("neb_ase_bfgs.traj", ":")
n_eval_ase = int(len(atoms_ase) - 2 * (len(atoms_ase)/n_images))
print("Number of function evaluations CI-NEB implemented in ASE:", n_eval_ase)

#print("Initial potential energy = ", neb_ase.get_potential_energy())
#print("Initial forces")
#print(neb_ase.get_forces())

for i,image in enumerate(neb_ase.images):
    print("%3d %18.10f" % (i, image.get_potential_energy()))
