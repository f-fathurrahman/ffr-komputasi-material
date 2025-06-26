import sys
sys.path.append('../')

from mini_mlip import MiniMLIP

#train_data = "ALL_ATOMS.xyz" # ALL_ATOMS_T500
#path_model = "LOGDIR_mini_Ni_fcc_352/" # need trailing /

#train_data = "ALL_ATOMS_T500.xyz"
#path_model = "LOGDIR_mini_Ni_fcc_T500/" # need trailing /

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

ff = MiniMLIP(descriptors=descriptor_dict, model=model_dict)
ff.run(mode="train", TrainData=train_data)
