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

#
#R_d_desc_alpha = desc.d_desc_dot_vec(R_d_desc, alphas.reshape(-1, dim_i))
#print("sum abs R_d_desc_alpha = ", np.sum(np.abs(R_d_desc_alpha)))
vecs = alphas.reshape(-1, dim_i)
print("vecs.shape now 1 = ", vecs.shape)
i, j = desc.tril_indices # these are atom indices
vecs = vecs.reshape(vecs.shape[0], -1, 3)
print("vecs.shape now 2 = ", vecs.shape)

print("R_d_desc.shape = ", R_d_desc.shape)
dvij = vecs[:,j,:] - vecs[:,i,:] # just to know the shape
print("dvij.shape = ", dvij.shape)

R_d_desc_alpha = np.einsum('...ij,...ij->...i', R_d_desc, vecs[:, j, :] - vecs[:, i, :])
print("R_d_desc_alpha.shape = ", R_d_desc_alpha.shape)
print("sum abs R_d_desc_alpha v2 = ", np.sum(np.abs(R_d_desc_alpha)))

"""
sum abs alphas = 227310078627.19122
sum abs R_d_desc_alpha =  70018251794.06943

vecs.shape now 1 =  (200, 27)
vecs.shape now 2 =  (200, 9, 3)
R_d_desc.shape =  (200, 36, 3)
dvij.shape =  (200, 36, 3)
R_d_desc_alpha.shape =  (200, 36)
sum abs R_d_desc_alpha v2 =  70018251794.06943

Console output

In [10]: np.dot(R_d_desc[0,0,:], dvij[0,0,:])
Out[10]: 2285994.9196487283

In [11]: R_d_desc_alpha[0,0]
Out[11]: 2285994.9196487283

In [12]: np.dot(R_d_desc[0,1,:], dvij[0,1,:])
Out[12]: -20410071.571386464

In [13]: R_d_desc_alpha[0,1]
Out[13]: -20410071.571386464

"""



#exit()


"""
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

R_train = task["R_train"]
E_pred, _ = gdml_predict.predict()
#for r in R_train:
#    print("r = ", r.shape)
#    E_pred, _ = gdml_predict.predict(r[None,:])
E_ref = np.squeeze(task['E_train'])

e_fact = np.linalg.lstsq(
    np.column_stack((E_pred, np.ones(E_ref.shape))), E_ref, rcond=-1
)[0][0]
corrcoef = np.corrcoef(E_ref, E_pred)[0, 1]

# Least squares estimate for integration constant.
c = np.sum(E_ref - E_pred) / E_ref.shape[0]
print("Recover integration constant: c = ", c)

#import matplotlib.pyplot as plt
#import matplotlib.style
#matplotlib.style.use("dark_background")
#
#sidx = np.argsort(E_ref)
#plt.plot(E_ref[sidx], label="E_ref")
#plt.plot(E_pred[sidx]+c, label="E_pred+c")
#plt.legend()
#plt.show()
"""

