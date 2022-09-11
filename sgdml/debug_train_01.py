import sys
import numpy as np
from my_sgdml.train import GDMLTrain

dataset = np.load("DATASET/ethanol_dft.npz")
n_train = 200

gdml_train = GDMLTrain()

"""
task = gdml_train.create_task(
    dataset, n_train, use_sym=False,
    valid_dataset=dataset, n_valid=1000,
    sig=20, lam=1e-10
)
"""

train_dataset = dataset
use_sym = False # set to False to use GDML not sGDML
valid_dataset = dataset
n_valid = 1000
sig = 20 # integer, kernel length scale
lam = 1e-10

# Default parameters
perms = None
use_E = True
use_E_cstr = False

use_E_cstr = use_E and use_E_cstr
print("use_E_cstr = ", use_E_cstr)

n_atoms = train_dataset["R"].shape[1]
print("n_atoms = ", n_atoms)


from my_sgdml.utils import io
md5_train = io.dataset_md5(train_dataset)
md5_valid = io.dataset_md5(valid_dataset)

print("md5_train = ", md5_train)
print("md5_valid = ", md5_valid)

from my_draw_strat_sample import my_draw_strat_sample

# if "E" in train_dataset:
print("Will call draw_strat_sample")
idxs_train = my_draw_strat_sample(train_dataset["E"], n_train)

excl_idxs = (
    idxs_train if md5_train == md5_valid else np.array([], dtype=np.uint)
)

idxs_valid = my_draw_strat_sample(
    valid_dataset["E"],
    n_valid,
    excl_idxs=excl_idxs,
)

R_train = train_dataset["R"][idxs_train, :, :]
task = {
    "type": "t",
    "code_version": "0.5.1",
    "dataset_name": train_dataset["name"].astype(str),
    "dataset_theory": train_dataset["theory"].astype(str),
    "z": train_dataset["z"],
    "R_train": R_train,
    "F_train": train_dataset["F"][idxs_train, :, :],
    "idxs_train": idxs_train,
    "md5_train": md5_train,
    "idxs_valid": idxs_valid,
    "md5_valid": md5_valid,
    "sig": sig,
    "lam": lam,
    "use_E": use_E,
    "use_E_cstr": use_E_cstr,
    "use_sym": use_sym,
}

if use_E:
    print("Using E_train in task")
    task["E_train"] = train_dataset["E"][idxs_train]

"""
# Not used
if "r_unit" in train_dataset and "e_unit" in train_dataset:
    print("Pass here 85")
    task["r_unit"] = train_dataset["r_unit"]
    task["e_unit"] = train_dataset["e_unit"]
"""

task['perms'] = np.arange(train_dataset['R'].shape[1])[None,:]  # no symmetries

#-----------------
# TRAIN START HERE
#-----------------


task = dict(task)  # make mutable

n_train, n_atoms = task['R_train'].shape[:2]

from my_sgdml.utils.desc import Desc
desc = Desc(n_atoms, max_processes=None)

n_perms = task['perms'].shape[0]
tril_perms = np.array([Desc.perm(p) for p in task['perms']])

dim_i = 3 * n_atoms
dim_d = desc.dim

perm_offsets = np.arange(n_perms)[:, None] * dim_d
tril_perms_lin = (tril_perms + perm_offsets).flatten('F')

lat_and_inv = None
callback = None

R = task['R_train'].reshape(n_train, -1)
R_desc, R_d_desc = desc.from_R(
    R,
    lat_and_inv=lat_and_inv,
    callback=partial(
        callback, disp_str='Generating descriptors and their Jacobians'
    )
    if callback is not None
    else None,
)


# Generate label vector.
E_train_mean = None
y = task['F_train'].ravel().copy()
#if task['use_E'] and task['use_E_cstr']:
#    print("Pass here in 134")
#    E_train = task['E_train'].ravel().copy()
#    E_train_mean = np.mean(E_train)
#    y = np.hstack((y, -E_train + E_train_mean))
#    # y = np.hstack((n*Ft, (1-n)*Et))

y_std = np.std(y)
print("y_std = ", y_std)
y /= y_std
