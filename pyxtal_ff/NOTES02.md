
Example good result (high R2), not validated:

```
from my_pyxtal_ff.pyxtal_ff_01 import PyXtal_FF
train_data = "DATASET_others/TiAl_gabung.xyz"
path_model = "LOGDIR_TiAl_gabung_v03/" # need trailing /

descriptor = {
    "Rc": 4.0, # make Rc smaller to speed up the calculation
    "type": "SO3",
    "parameters": {"nmax": 4, "lmax": 3},
    "ncpu": 1,
}

model = {
    "system" : ["Ti", "Al"],
    "hiddenlayers": [30, 30],
    "force_coefficient": 0.01, # Need to tune this
    "epoch": 100,
    "batch_size": 32,
    "path": path_model,
    "memory": "out",
    "optimizer": {"method": "LBFGS"}
}

ff = PyXtal_FF(descriptors=descriptor, model=model)
ff.run(mode="train", TrainData=train_data)
```

It seems that stress fitting is not supported.

