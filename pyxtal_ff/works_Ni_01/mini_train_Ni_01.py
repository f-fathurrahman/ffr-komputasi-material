import sys
sys.path.append('../')

from mini_mlip import MiniMLIP

train_data = "ALL_ATOMS.xyz" # ALL_ATOMS_T500
path_model = "LOGDIR_mini_Ni_fcc/" # need trailing /

descriptor = {
    "Rc": 4.0,
    "type": "SO3",
    "parameters": {
        "nmax": 4,
        "lmax": 3
    },
    "ncpu": 1,
}

model = {
    "algorithm": "PR",
    "system" : ["Ni"],
    "path": path_model,
    'force_coefficient': 0.001,
}

ff = MiniMLIP(descriptors=descriptor, model=model)
ff.run(mode="train", TrainData=train_data)
