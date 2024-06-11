import numpy as np
import scipy
import os

import numpy.random
numpy.random.seed(1234)

from prepare_data import import_from_pickle

R_desc, R_d_desc, tril_perms_lin, desc, task, y = import_from_pickle()

sig = 20
lam = 1e-10

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
# sum abs K = 7919.284079790718
# sum abs y = 4001.788905789174


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
# 12315.115002408249 WSL2 AMD Ryzen
# 12315.002994660537 Ubuntu 20.04 native intel i7
# 12315.274852764584 Windows 11

print(f"sum abs alphas = {np.sum(np.abs(alphas))}")



y = task["F_train"].ravel().copy()
y_std = np.std(y)

dim_i = desc.dim_i

R_d_desc_alpha = desc.d_desc_dot_vec(R_d_desc, alphas.reshape(-1, dim_i))
print("sum abs R_d_desc_alpha = ", np.sum(np.abs(R_d_desc_alpha)))

# Prepare model dict
model = {
    'type': 'm',
    'code_version': '0.0.1ffr',
    'dataset_name': task['dataset_name'],
    'dataset_theory': task['dataset_theory'],
    'solver_name': 'analytic',
    'z': task['z'],
    'idxs_train': task['idxs_train'],
    'idxs_valid': task['idxs_valid'],
    'n_test': 0,
    'md5_test': None,
    'f_err': {'mae': np.nan, 'rmse': np.nan},
    'c': 0.0, # not yet calculated
    'std': y_std,
    'sig': task['sig'],
    'lam': task['lam'],
    'alphas_F': alphas,
    'R_desc': R_desc.T,   # need transpose!!!!
    'R_d_desc_alpha': R_d_desc_alpha,
    'perms': task['perms'],
    'tril_perms_lin': tril_perms_lin,
    'use_E': task['use_E'],
}


from my_predict import GDMLPredict
gdml_predict = GDMLPredict(
    model,
    max_processes=1,
    log_level=0,
)

gdml_predict.set_R_desc(R_desc)
gdml_predict.set_R_d_desc(R_d_desc)

# Predict, one data point

R_train = task["R_train"]
r = R_train[9]
print("r original shape = ", r.shape)
Natoms = r.shape[0]
print("Natoms = ", Natoms)
r = r.reshape(1,Natoms*3) # we need to reshape it

if False:
    E_pred, F_pred = gdml_predict.predict(r)
    F_pred = F_pred.reshape(Natoms,3)
    print("E_pred = ", E_pred)
    print("F_pred = ")
    print(F_pred)
    exit()

# Compute descriptor and its Jacobian
lat_and_inv = None
r_desc, r_d_desc = desc.from_R(r, lat_and_inv)

print("sum abs r_desc = ", np.sum(np.abs(r_desc)))


# These are only for the case of n_perms == 1


diff_ab = np.subtract(
    np.broadcast_to(r_desc, R_desc.shape),
    R_desc
)
print("sum abs diff_ab = ", np.sum(np.abs(diff_ab)))

diag_scale_fact = 5.0 / sig

norm_ab = np.sqrt(5) * np.linalg.norm(diff_ab, axis=1)
print("sum abs norm_ab = ", np.sum(np.abs(norm_ab)))

mat52_base_fact = 5.0 / (3 * sig ** 3)
mat52_base = np.exp(-norm_ab/sig)
mat52_base *= mat52_base_fact
print("sum abs mat52_base = ", np.sum(np.abs(mat52_base)))

# column wise dot product
a_x2 = np.einsum("ji,ji->j", diff_ab, R_d_desc_alpha)
print("sum abs a_x2 = ", np.sum(np.abs(a_x2)))

F = (a_x2 * mat52_base).dot(diff_ab) * diag_scale_fact
print("sum abs F before after dot = ", np.sum(np.abs(F)))

mat52_base *= norm_ab + sig
print("sum abs mat52_base after scaling ", np.sum(np.abs(mat52_base)))

F -= mat52_base.dot(R_d_desc_alpha)
print("sum abs F before after dot with R_d_desc_alpha = ", np.sum(np.abs(F)))

E_pred0 = a_x2.dot(mat52_base)*y_std
print("E_pred0 = ", E_pred0)

# Multiply with r_d_desc with F
#r_d_F = desc.vec_dot_d_desc(r_d_desc, F)*y_std

r_jac = np.copy(r_d_desc)
vecs = np.copy(F)

if r_jac.ndim == 2:
    r_jac = r_jac[None, ...]

if vecs.ndim == 1:
    vecs = vecs[None, ...]

assert (
    r_jac.shape[0] == 1
    or vecs.shape[0] == 1
    or R_d_desc.shape[0] == vecs.shape[0]
)
# either multiple descriptors or multiple vectors at once, not both
# (or the same number of both, than it will must be a multidot) (???)

n = np.max((r_jac.shape[0], vecs.shape[0]))
# In the present case: n = 1
i, j = desc.tril_indices


"""
In [13]: vecs[...,None].shape
Out[13]: (1, 36, 1)

In [14]: vecs.shape
Out[14]: (1, 36)
"""

out = np.zeros((n, Natoms, Natoms, 3))
out[:, i, j, :] = r_jac * vecs[..., None]
out[:, j, i, :] = -out[:, i, j, :]
r_d_F = out.sum(axis=1).reshape(n, -1) * y_std

F_pred = r_d_F.reshape(Natoms,3)
print(F_pred)

