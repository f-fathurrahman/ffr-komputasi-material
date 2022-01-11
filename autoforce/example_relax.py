# +
import numpy as np
from ase.atoms import Atoms
from ase.calculators.vasp import Vasp
from ase.calculators.emt import EMT
from ase.optimize import LBFGS, BFGS
from util_flake import generate_random_cluster
from calculator_active import ActiveCalculator, kcal_mol

import torch
torch.set_default_tensor_type(torch.DoubleTensor)

np.random.seed(1234)

# random cluster generation
#ngold = 20
#min_dist = 2.5
#positions = generate_random_cluster(ngold, min_dist)
#atoms = Atoms(numbers=ngold*[79], positions=positions)
#atoms.center(vacuum=5.)
#atoms.pbc = True
#atoms.write('INITIAL.xyz', format='extxyz')

import ase.io
atoms = ase.io.read("INITIAL.xyz")
print(atoms.pbc)
print("Atoms read successfully")
#print(atoms.positions)
#exit()


# Ab initio calculator; for now we just use EMT instead of vasp
# abinitio = Vasp(command="mpirun -n 16 vasp_std", directory='vasp')
abinitio = EMT()


# ML calculator
active_kwargs = {
    'calculator': abinitio,
    'ediff': 0.1*kcal_mol,  # decrease for more accuracy but lower speed
    'fdiff': 0.1*kcal_mol,  # decrease for more accuracy but lower speed
     'kernel_kw': {'cutoff': 6.5, 'lmax': 3, 'nmax': 3},
    # 'kernel_kw': {'cutoff': 6., 'lmax': 3, 'nmax': 3, 'species': [79]}, # <- faster
    # 'veto': {'forces': 8.}  # for vetoing ML updates for very high energy structures
}
calc = ActiveCalculator(**active_kwargs)
atoms.calc = calc

#atoms.calc = abinitio # Using EAM


# relax
maxforce = 0.01
dyn = LBFGS(atoms, trajectory='relax.traj')
#dyn = BFGS(atoms, trajectory='relax.traj')
#dyn.run(fmax=maxforce, steps=10)
dyn.run(fmax=maxforce)

#exit()

atoms.write('optimized.xyz', format='extxyz')


# For history dependent algorithms such as LBFGS, if ML
# updates becomes problematic, dyn.run can be replaced by:
#
# for _ in dyn.irun(fmax=maxforce):
#    if calc.updated:
#        dyn.initialize()    # clears history


# Getting closer to the local minima by forceful updating
# the model; this will cause at least 2 more ab initio
# calculations (one of which is needed anyway to determine
# the proximity to the actual minima).
while True:
    if calc.update_data(try_fake=False):
        calc.update(data=False)
        calc.results.clear()
        dyn.initialize()
        dyn.run(fmax=maxforce)
    else:
        break


# Calculate exact energy & forces for the final coordinates.
# The optimized coordinates and ab initio energy & forces
# will be written to "active_FP.traj".
energy, forces = atoms.calc._test()
f_rms = np.sqrt(np.mean(forces**2))
f_max = abs(forces).max()
report = f"""
    relaxation result:
    energy:      {energy}
    force (rms): {f_rms}
    force (max): {f_max}
"""
print(report)


# save the final structure
# atoms.write('CONTCAR')
atoms.write('optimized.xyz', format='extxyz')
