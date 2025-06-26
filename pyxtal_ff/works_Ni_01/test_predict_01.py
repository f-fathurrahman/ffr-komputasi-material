import sys
sys.path.append('../')

from mini_mlip import MiniMLIP, MiniMLIPCalculator, mini_mlip_optimize

#train_data = "ALL_ATOMS.xyz"
#path_model = "LOGDIR_mini_Ni_fcc_352/" # need trailing /

train_data = "DATASET_26062025/DATA_Ni_26_06_2025.xyz"
path_model = "LOGDIR_mini_26062025/" # need trailing /


descriptor_dict = {
    "Rc": 4.0,
    "type": "SO3",
    "parameters": {
        "nmax": 4,
        "lmax": 3
    },
    "ncpu": 1,
    "stress": False,
}

model_dict = {
    "algorithm": "PR",
    "system" : ["Ni"],
    "path": path_model,
    "force_coefficient": 0.001,
    "stress_coefficient": None,
}

ff = MiniMLIP(model=model_dict, descriptors=descriptor_dict)
ff.run(mode='predict', mliap=path_model+"PolyReg-checkpoint.pth")
calc = MiniMLIPCalculator(ff=ff)

from ase.build import bulk
from ase import units
atoms = bulk('Ni', 'fcc', cubic=True)
atoms.calc = calc
print('\ninitial cell para: ', atoms.get_cell())
print('initial energy: ', atoms.get_potential_energy())
print('initial stress', -atoms.get_stress()/units.GPa)

v1 = atoms.get_cell()[:,0]
import numpy as np
print("v1.norm = ", np.linalg.norm(v1))

# geometry optimization
atoms = mini_mlip_optimize(atoms, box=True)
print('\nequlibrium cell para: ', atoms.get_cell())
print('equlirium energy: ', atoms.get_potential_energy())
print('equlibrium stress', -atoms.get_stress()/units.GPa)
v1 = atoms.get_cell()[:,0]
import numpy as np
print("v1.norm = ", np.linalg.norm(v1))