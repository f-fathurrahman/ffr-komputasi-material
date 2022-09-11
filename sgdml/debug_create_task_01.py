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
# convert to matrix
