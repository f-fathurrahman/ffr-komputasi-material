# Config parameters
precision       = "default"       # "default"(fp32), "low"(mixed 32-16), "high"(fp64)
save_name       = "model.pkl"     # model save path
model_type      = "energy"        # "energy" or "atomic" (e.g. wannier)
atomic_sel      = [0]             # select atom type for prediction (only for "atomic" model)
atomic_label    = "atomic_dipole" # data file prefix for "atomic" model; string must contain "atomic"

# Dataset in DeepMD-kit format; nested paths like [[dat1,dat2],[dat3]] allowed
# Note: Here the atomic type index of dat1,dat2 must be the same, but that of dat3 can be different
# train_paths     = ["data/chunyi_dplr/data/dipole_data"]
train_paths     = ["DATASET/CH4_data/training_data"]  # paths to training data
use_val_data    = False           # if not, next line is ignored
val_paths       = ["DATASET/CH4_data/validata_data"]    # paths to validation data

# Model parameters
rcut            = 6.0             # cutoff radius (Angstrom)
use_2nd_tensor  = True            # Use 2nd order tensor descriptor for more accuracy, slightly slower
use_mp          = False           # Use message passing (DP-MP) model for even more accuracy (slower) 
compress        = True            # Compress model after training for faster inference. Rec: True  
embed_widths    = [32,32,64]      # Rec: [32,32,64] (for accuracy try [48,48,96])
embedMP_widths  = [32,64]         # Rec: [32,64]; Only used in MP; (Try [64,64] or [96,96] according to embed_widths)
fit_widths      = [64,64,64]      # For "atomic" model, fit_widths[-1] must equal embed_widths[-1](DP)/embedMP_widths[-1](DP-MP)
axis_neurons    = 12              # Rec: 8-16

# Training parameters
label_bs        = 160             # training batch size in Nlabels, Rec: 128-512; Nframes_per_batch = ceil(label_bs/Nlabels_per_frame)
val_label_bs    = 1024            # validation batch size in Nlabels. Too much can cause OOM error.
lr              = 0.001           # learning rate at start. Rec: 0.001/0.002 for "energy", 0.01 for "atomic"
s_pref_e        = 0.02            # starting prefactor for energy loss
l_pref_e        = 1               # limit prefactor for energy loss, increase for energy accuracy
s_pref_f        = 1000            # starting prefactor for force loss
l_pref_f        = 1               # limit prefactor for force loss, increase for force accuracy
total_steps     = 30000           # total training steps. Rec: 1e6 for "energy", 1e5 for "atomic"
print_every     = 1000            # for printing loss and validation

# parameters you usually don"t need to change
lr_limit        = 5e-7            # learning rate at end of training
compress_Ngrids = 1024            # Number of intervals used in compression
compress_rmin   = 0.6             # Lower bound for interatomic distance in compression
beta2           = 0.99            # adam optimizer parameter
l_smoothing     = 20              # smoothing factor for loss printing
decay_steps     = 5000            # learning rate exponentially decays every decay_steps
getstat_bs      = 64              # batch size for computing model statistics at initialization

# From here on you don"t need to change anything unless you know what you are doing
import jax
import datetime
from time import time

from my_deepmd_jax import data

# This is not yet needed for loading data
if precision == "default":
    print("Using default precision float32")
    jax.config.update("jax_default_matmul_precision", "float32")
if precision == "high":
    print("Using float64")
    jax.config.update("jax_enable_x64", True)


TIC = time()
print("# Program start at", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "on device:", jax.devices()[:1])

if model_type == "energy":
    labels = ["coord","box"] + ["force","energy"]
else:
    labels = ["coord","box"] + [atomic_label]

print("labels = ", labels)

# Load data
train_data = data.DPDataset(train_paths, labels, {"atomic_sel": atomic_sel})

if use_val_data:
    val_data = data.DPDataset(val_paths, labels, {"atomic_sel": atomic_sel})
else:
    val_data = None


import numpy as np
from jax import vmap
import jax.numpy as jnp

boxes = train_data.data["box"]

print("boxes = ", boxes)

# def compute_lattice_candidate(boxes, rcut): # boxes (nframes,3,3)
N = 2  # This algorithm is heuristic and subject to change. Increase N in case of missing neighbors.

ortho = not vmap(lambda box: box - jnp.diag(jnp.diag(box)))(boxes).any()
# ortho is a scalar?

recp_norm = jnp.linalg.norm((jnp.linalg.inv(boxes)), axis=-1)    # (nframes,3)
n = np.ceil(rcut * recp_norm - 0.5).astype(int).max(0)           # (3,)
lattice_cand = jnp.stack(
    np.meshgrid(
        range(-n[0],n[0]+1),
        range(-n[1],n[1]+1),
        range(-n[2],n[2]+1),
        indexing='ij'
    ),
    axis=-1
).reshape(-1,3)

trial_points = jnp.stack(
    np.meshgrid(
        np.arange(-N,N+1),
        np.arange(-N,N+1),
        np.arange(-N,N+1)
    ),
    axis=-1
).reshape(-1,3) / (2*N)

is_neighbor = jnp.linalg.norm(
    (lattice_cand[:,None] - trial_points)[None] @ boxes[:,None],
axis=-1) < rcut  # (nframes,l,t)

lattice_cand = np.array(lattice_cand[is_neighbor.any((0,2))])

lattice_max = is_neighbor.sum(1).max().item()


print('Lattice vectors for neighbor images: Max %d out of %d condidates.' % (lattice_max, len(lattice_cand)))
dict_res = {
    'lattice_cand': tuple(map(tuple, lattice_cand)),
    'lattice_max': lattice_max,
    'ortho': ortho
}