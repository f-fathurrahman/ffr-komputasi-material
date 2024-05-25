from functools import partial
import os

import numpy as np
import scipy as sp
import timeit
from gdml_train_function import my_assemble_kernel_mat

def my_analytic_solve(desc, task, R_desc, R_d_desc, tril_perms_lin, y):

    sig = task['sig']
    lam = task['lam']
    use_E_cstr = task['use_E_cstr']

    n_train, dim_d = R_d_desc.shape[:2]
    n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)
    dim_i = 3 * n_atoms

    K = -my_assemble_kernel_mat(
        R_desc,
        R_d_desc,
        tril_perms_lin,
        sig,
        desc,
        use_E_cstr=use_E_cstr
    )  # Flip sign to make convex
    print("shape of kernel matrix K = ", K.shape)

    start = timeit.default_timer()

    if K.shape[0] == K.shape[1]:

        K[np.diag_indices_from(K)] += lam  # Regularize

        print('Solving linear system (Cholesky factorization)')

        try:

            # Cholesky (do not overwrite K in case we need to retry)
            L, lower = sp.linalg.cho_factor(
                K, overwrite_a=False, check_finite=False
            )
            alphas = -sp.linalg.cho_solve(
                (L, lower), y, overwrite_b=False, check_finite=False
            )

        except np.linalg.LinAlgError:  # Try a solver that makes less assumptions
            #
            print('Solving linear system (LU factorization)')
            try:
                # LU
                alphas = -sp.linalg.solve(
                    K, y, overwrite_a=True, overwrite_b=True, check_finite=False
                )
            except MemoryError:
                print('Not enough memory to train this system using a closed form solver.')
                print()
                os._exit(1)

        except MemoryError:
            print('Not enough memory to train this system using a closed form solver.')
            print()
            os._exit(1)
    else:

        print('Solving over-determined linear system (least squares approximation)')
        # Least squares for non-square K
        alphas = -np.linalg.lstsq(K, y, rcond=-1)[0]

    stop = timeit.default_timer()

    print(f"Finished training in {stop-start} s")

    return alphas

