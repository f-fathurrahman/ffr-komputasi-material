import INPUT01 as INPUT

# From here on you don"t need to change anything unless you know what you are doing
import jax
import datetime
from time import time

from my_deepmd_jax import data

# This is not yet needed for loading data
if INPUT.precision == "default":
    print("Using default precision float32")
    jax.config.update("jax_default_matmul_precision", "float32")
if INPUT.precision == "high":
    print("Using float64")
    jax.config.update("jax_enable_x64", True)


TIC = time()
print("# Program start at", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "on device:", jax.devices()[:1])

if INPUT.model_type == "energy":
    labels = ["coord","box"] + ["force","energy"]
else:
    labels = ["coord","box"] + [INPUT.atomic_label]

print("labels = ", labels)

# Load data
train_data = data.DPDataset(INPUT.train_paths, labels, {"atomic_sel": INPUT.atomic_sel})

if INPUT.use_val_data:
    val_data = data.DPDataset(INPUT.val_paths, labels, {"atomic_sel": INPUT.atomic_sel})
else:
    val_data = None

# Compute lattice vectors?
train_data.compute_lattice_candidate(INPUT.rcut)
if INPUT.use_val_data:
    val_data.compute_lattice_candidate(INPUT.rcut)
