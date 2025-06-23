import sys
sys.path.append('../')

from mini_mlip import MiniMLIP, MiniMLIPCalculator

train_data = "ALL_ATOMS.xyz"
path_model = "LOGDIR_mini_Ni_fcc/" # need trailing /

descriptor_dict = {
    "Rc": 4.0,
    "type": "SO3",
    "parameters": {
        "nmax": 4,
        "lmax": 3
    },
    "ncpu": 1,
}

model_dict = {
    "algorithm": "PR",
    "system" : ["Ni"],
    "path": path_model,
    'force_coefficient': 0.001,
}

ff = MiniMLIP(model=model_dict, descriptors=descriptor_dict)
ff.run(mode='predict', mliap=path_model+"PolyReg-checkpoint.pth")
calc = MiniMLIPCalculator(ff=ff)