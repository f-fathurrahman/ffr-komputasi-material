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
n = np.ceil(INPUT.rcut * recp_norm - 0.5).astype(int).max(0)           # (3,)
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
axis=-1) < INPUT.rcut  # (nframes,l,t)

lattice_cand = np.array(lattice_cand[is_neighbor.any((0,2))])

lattice_max = is_neighbor.sum(1).max().item()

print('Lattice vectors for neighbor images: Max %d out of %d condidates.' % (lattice_max, len(lattice_cand)))
dict_res = {
    'lattice_cand': tuple(map(tuple, lattice_cand)),
    'lattice_max': lattice_max,
    'ortho': ortho
}
