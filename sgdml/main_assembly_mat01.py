import numpy as np
import scipy
import os

import numpy.random
numpy.random.seed(1234)

from prepare_data import import_from_pickle

R_desc, R_d_desc, tril_perms_lin, desc, task, y = import_from_pickle()

sig = 20
lam = 1e-10

# Run this to compare with the original version which uses multiprocessing
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

print(f"sum abs K = {np.sum(np.abs(K))}")
print(f"sum abs y = {np.sum(np.abs(y))}")

#
# This is the analytic.solve step
#
if K.shape[0] == K.shape[1]:
    K[np.diag_indices_from(K)] += lam  # Regularize
    try:
        # To force using LU decomposition
        #raise np.linalg.LinAlgError
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

print("average alphas = ", np.average(alphas))
# 12315.115002408249

print(f"sum abs alphas = {np.sum(np.abs(alphas))}")