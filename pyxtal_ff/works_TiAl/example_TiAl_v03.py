from pyxtal_ff_01 import PyXtal_FF

train_data = "DATASET_others/TiAl_gabung.xyz"
path_model = "LOGDIR_TiAl_gabung_acsf/" # need trailing /

symmetry = {
    "G2": {
        "eta": [0.035709, 0.071418, 0.178545, 0.35709, 0.71418, 1.78545],
        "Rs": [0]
    },
    "G4": {
        "lambda": [-1, 1],
        "zeta": [1],
        "eta": [0.035709, 0.071418, 0.178545, 0.35709]
    }
}

descriptor = {
    "type": "BehlerParrinello",
    "parameters": symmetry,
    "Rc": 5.0,
    "ncpu": 1,
}

model = {
    "system" : ["Ti", "Al"],
    "hiddenlayers": [30, 30],
    "force_coefficient": 0.0001,
    "epoch": 50,
    "batch_size": 32,
    "path": path_model,
    "memory": "out",
    "optimizer": {"method": "LBFGS"}
}

ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode="train", TrainData=train_data)
