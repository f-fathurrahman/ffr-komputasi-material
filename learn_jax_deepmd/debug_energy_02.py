import INPUT01 as INPUT

#from INPUT01 import *

# override some input parameters
INPUT.total_steps = 1
INPUT.print_every = 1


# From here on you don"t need to change anything unless you know what you are doing
import numpy as np
from jax import jit, random, tree_util
import jax, optax, datetime
import flax.linen as nn
from time import time
from functools import partial

from my_deepmd_jax import data, utils
from my_deepmd_jax.dpmodel import DPModel


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
print("Program start at", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "on device:", jax.devices()[:1])

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


params = {
    "embed_widths": INPUT.embed_widths[:-1] if INPUT.use_mp else INPUT.embed_widths,
    "embedMP_widths": INPUT.embed_widths[-1:] + INPUT.embedMP_widths if INPUT.use_mp else None,
    "fit_widths": INPUT.fit_widths,
    "axis": INPUT.axis_neurons,
    "Ebias": train_data.fit_energy() if INPUT.model_type == "energy" else None,
    "rcut": INPUT.rcut,
    "use_2nd": INPUT.use_2nd_tensor,
    "use_mp": INPUT.use_mp,
    "atomic":True if INPUT.model_type == "atomic" else False,
    "nsel": INPUT.atomic_sel if INPUT.model_type == "atomic" else None,
    "out_norm": 1. if INPUT.model_type == "energy" else train_data.get_atomic_label_scale()
}

data_stats = train_data.get_stats(INPUT.rcut, INPUT.getstat_bs)
model = DPModel(params|data_stats) # | will concat two dicts
print("Model params:", model.params)


# initialize model variables
batch, type_count, lattice_args = train_data.get_batch(1) # get one data point
static_args  = nn.FrozenDict({"type_count": type_count, "lattice": lattice_args})
variables = model.init(
    random.PRNGKey(np.random.randint(42)),
    batch["coord"][0],
    batch["box"][0],
    static_args
)


coord_N3 = batch["coord"][0]
box_33 = batch["box"][0]
nbrs_nm = None # not used?

model.apply({"params": variables["params"]}, coord_N3, box_33, static_args)


