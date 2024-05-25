import numpy as np
import scipy
import os

import numpy.random
numpy.random.seed(1234)

from my_sgdml.train import GDMLTrain

filename = "DATASET/ethanol_dft.npz"

dataset = np.load(filename)
n_train = 200
gdml_train = GDMLTrain()

#
# We create task here (manually)
#
train_dataset = dataset
use_sym = False # set to False to use GDML not sGDML
valid_dataset = dataset
n_valid = 1000 # not used
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

#
# Sampling the dataset
#
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
F_train = train_dataset["F"][idxs_train, :, :]
task = {
    "type": "t",
    "code_version": "0.5.1",
    "dataset_name": train_dataset["name"].astype(str),
    "dataset_theory": train_dataset["theory"].astype(str),
    "z": train_dataset["z"],
    "R_train": R_train,
    "F_train": F_train,
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

task["perms"] = np.arange(train_dataset["R"].shape[1])[None,:]  # no symmetries

# Now task is ready



#-----------------
# TRAIN START HERE
#-----------------

task = dict(task)  # make mutable (need this?)

n_train, n_atoms = task["R_train"].shape[:2]

from my_sgdml.utils.desc import Desc
desc = Desc(n_atoms, max_processes=None)

n_perms = task["perms"].shape[0]
tril_perms = np.array([Desc.perm(p) for p in task["perms"]])

dim_i = 3 * n_atoms
dim_d = desc.dim

perm_offsets = np.arange(n_perms)[:, None] * dim_d
tril_perms_lin = (tril_perms + perm_offsets).flatten("F")

lat_and_inv = None
callback = None

R = task["R_train"].reshape(n_train, -1)
R_desc, R_d_desc = desc.from_R(
    R,
    lat_and_inv=lat_and_inv,
    callback=None,
)


# Generate label vector.
E_train_mean = None
y = task["F_train"].ravel().copy()
y_std = np.std(y)
print("y_std = ", y_std)
y /= y_std


"""
from my_sgdml.solvers.analytic import Analytic
analytic = Analytic(gdml_train, desc, callback=callback)
alphas = analytic.solve(task, R_desc, R_d_desc, tril_perms_lin, y)
print(alphas.shape)
"""

sig = task["sig"]
lam = task["lam"]
use_E_cstr = task["use_E_cstr"]

n_train, dim_d = R_d_desc.shape[:2]
n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)
dim_i = 3 * n_atoms

print("Assembling kernel matrix")

K = -gdml_train._assemble_kernel_mat(
    R_desc,
    R_d_desc,
    tril_perms_lin,
    sig,
    desc,
    use_E_cstr=False,
    callback=None,
)  # Flip sign to make convex

print(K.shape)

#
# This is the analytic.solve step
#
if K.shape[0] == K.shape[1]:
    K[np.diag_indices_from(K)] += lam  # Regularize
    try:
        print("Trying Cholesky decomposition")
        # Cholesky (do not overwrite K in case we need to retry)
        L, lower = scipy.linalg.cho_factor(
            K, overwrite_a=False, check_finite=False
        )
        alphas = -scipy.linalg.cho_solve(
            (L, lower), y, overwrite_b=False, check_finite=False
        )
        print("Solving linear equation successful")
    except np.linalg.LinAlgError:  # Try a solver that makes less assumptions
        try:
            print("Trying LU decomposition")
            # LU
            alphas = -scipy.linalg.solve(
                K, y, overwrite_a=True, overwrite_b=True, check_finite=False
            )
        except MemoryError:
            print("Not enough memory to train this system using a closed form solver.")
            print()
            os._exit(1)
    except MemoryError:
        print("Not enough memory to train this system using a closed form solver.")
        print()
        os._exit(1)
else:
    print("Using least squares")
    # Least squares for non-square K
    alphas = -np.linalg.lstsq(K, y, rcond=-1)[0]



# The parameters are found, now we need to create a model
# and find

alphas_E = None
alphas_F = alphas

solver_keys = {} # remove this
print("After finding the parameters, create model")
print("solver_keys = ", solver_keys)
model = gdml_train.create_model(
    task,
    "analytic",
    R_desc,
    R_d_desc,
    tril_perms_lin,
    y_std,
    alphas_F,
    alphas_E=alphas_E,
)
model.update(solver_keys)

# Recover integration constant.
# Note: if energy constraints are included in the kernel (via "use_E_cstr"), do not
# compute the integration constant, but simply set it to the mean of the training energies
# (which was subtracted from the labels before training).
if model["use_E"]:
    c = (
        gdml_train._recov_int_const(model, task, R_desc=R_desc, R_d_desc=R_d_desc)
        if E_train_mean is None
        else E_train_mean
    )
    print("Recover integration constant: c = ", c)
    if c is None:
        # Something does not seem right.
        # Turn off energy predictions for this model, only output force predictions.
        model["use_E"] = False
    else:
        model["c"] = c

