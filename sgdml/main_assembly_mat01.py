import numpy as np
import scipy
import os

import numpy.random
numpy.random.seed(1234)

from my_draw_strat_sample import my_draw_strat_sample
from my_desc_from_R import my_desc_from_R
from my_sgdml.utils import io
from my_sgdml.utils.desc import Desc


# Prepare input for assembly_kernel_mat
def prepare_data(filename="DATASET/ethanol_dft.npz", n_train=200):

    dataset = np.load(filename)

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

    md5_train = io.dataset_md5(train_dataset)
    md5_valid = io.dataset_md5(valid_dataset)

    print("md5_train = ", md5_train)
    print("md5_valid = ", md5_valid)

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

    # from training step

    task = dict(task)  # make mutable (need this?)
    n_train, n_atoms = task["R_train"].shape[:2]

    desc = Desc(n_atoms, max_processes=1)

    n_perms = task["perms"].shape[0]
    tril_perms = np.array([Desc.perm(p) for p in task["perms"]])

    dim_d = desc.dim

    perm_offsets = np.arange(n_perms)[:, None] * dim_d
    tril_perms_lin = (tril_perms + perm_offsets).flatten("F")

    lat_and_inv = None
    R = task['R_train'].reshape(n_train, -1)
    R_desc, R_d_desc = my_desc_from_R(
        desc, R, lat_and_inv=lat_and_inv
    )



    # Generate label vector.
    y = task["F_train"].ravel().copy()
    y_std = np.std(y)
    y /= y_std

    sig = task["sig"]
    lam = task["lam"]
    use_E_cstr = task["use_E_cstr"]

    n_train_, dim_d = R_d_desc.shape[:2]
    assert n_train_ == n_train
    n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)

    return R_desc, R_d_desc, tril_perms_lin, desc, task, y



R_desc, R_d_desc, tril_perms_lin, desc, task, y = prepare_data()

sig = 20
lam = 1e-10

#from gdml_train_function_02 import my_assemble_kernel_mat
#K = -my_assemble_kernel_mat(
#    R_desc,
#    R_d_desc,
#    tril_perms_lin,
#    sig,
#    desc,
#    use_E_cstr=False
#)  # Flip sign to make convex
#print("shape of kernel matrix K = ", K.shape)
#np.savez("ORIG_K_ethanol.npz", K)

n_train, dim_d = R_d_desc.shape[:2]
dim_i = 3 * int((1 + np.sqrt(8 * dim_d + 1)) / 2)

# Determine size of kernel matrix.
K_n_rows = n_train * dim_i
K_n_cols = K_n_rows


K = np.zeros((K_n_rows,K_n_cols))
from debug_assemble_kernel_mat_wkr import debug_assemble_kernel_mat_wkr
for j in range(n_train):
    debug_assemble_kernel_mat_wkr(
        K, R_desc, R_d_desc, desc, j, tril_perms_lin, sig )
K *= -1 # flip sign



"""
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

print("average alpha = ", np.average(alphas))
# 12315.115002408249
"""