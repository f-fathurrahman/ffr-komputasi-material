import multiprocessing as mp
Pool = mp.get_context('fork').Pool

import timeit
import numpy as np

from functools import partial

# Global variable, will be shared among functions within this file only (?)
glob = {}


# HUH?
def _share_array(arr_np, typecode_or_type):
    arr = mp.RawArray(typecode_or_type, arr_np.ravel())
    return arr, arr_np.shape




def my_assemble_kernel_mat(
    R_desc,
    R_d_desc,
    tril_perms_lin,
    sig,
    desc,  # TODO: document me
    use_E_cstr=False,
    col_idxs=np.s_[:],  # TODO: document me
    alloc_extra_rows=0,  # TODO: document me
):

    global glob

    # Note: This function does not support unsorted (ascending) index arrays.
    # if not isinstance(col_idxs, slice):
    #    assert np.array_equal(col_idxs, np.sort(col_idxs))

    n_train, dim_d = R_d_desc.shape[:2]
    dim_i = 3 * int((1 + np.sqrt(8 * dim_d + 1)) / 2)

    # Determine size of kernel matrix.
    K_n_rows = n_train * dim_i
    K_n_cols = K_n_rows

    print("K_n_rows = ", K_n_rows)
    print("K_n_cols = ", K_n_cols)

    # Account for additional rows and columns due to energy constraints in the kernel matrix.
    if use_E_cstr:
        K_n_rows += n_train
        K_n_cols += n_train

    # Make sure no indices are outside of the valid range.
    if K_n_cols > K_n_rows:
        raise ValueError('Columns indexed beyond range.')

    print("K_n_rows = ", K_n_rows, " K_n_cols = ", K_n_cols)
    print("alloc_extra_rows = ", alloc_extra_rows)

    # HUH???
    K = mp.RawArray('d', (K_n_rows + alloc_extra_rows) * K_n_cols)
    glob['K'], glob['K_shape'] = K, (K_n_rows + alloc_extra_rows, K_n_cols)
    glob['R_desc'], glob['R_desc_shape'] = _share_array(R_desc, 'd')
    glob['R_d_desc'], glob['R_d_desc_shape'] = _share_array(R_d_desc, 'd')

    glob['desc_func'] = desc

    print("type(K) = ", type(K))

    start = timeit.default_timer()

    for j in range(n_train):
        my_assemble_kernel_mat_wkr(
            j, tril_perms_lin, sig, use_E_cstr=False,
            exploit_sym=True, cols_m_limit=None
        )

    stop = timeit.default_timer()

    print(f"Time to assemble: {stop - start}")

    # Release some memory.
    glob.pop('K', None)
    glob.pop('R_desc', None)
    glob.pop('R_d_desc', None)

    # return np.frombuffer(K).reshape(glob['K_shape'])
    return np.frombuffer(K).reshape((K_n_rows + alloc_extra_rows), K_n_cols)




def my_assemble_kernel_mat_wkr(
    j, tril_perms_lin, sig, use_E_cstr=False, exploit_sym=False, cols_m_limit=None
):

    #print("*** ENTER _assemble_kernel_mat_wkr for index j = ", j)

    global glob

    # XXX R_desc, R_d_desc, and K are communicated via global dictionary glob
    R_desc = np.frombuffer(glob['R_desc']).reshape(glob['R_desc_shape'])
    R_d_desc = np.frombuffer(glob['R_d_desc']).reshape(glob['R_d_desc_shape'])
    K = np.frombuffer(glob['K']).reshape(glob['K_shape'])
    # XXX
    desc_func = glob['desc_func']

    n_train, dim_d = R_d_desc.shape[:2]
    n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)
    dim_i = 3 * n_atoms
    n_perms = int(len(tril_perms_lin) / dim_d)

    blk_j = slice(j * dim_i, (j + 1) * dim_i)
    keep_idxs_3n = slice(None)  # same as [:]

    #print("blk_j = ", blk_j)

    # TODO: document this exception
    if use_E_cstr and not (cols_m_limit is None or cols_m_limit == n_train):
        raise ValueError(
            '\'use_E_cstr\'- and \'cols_m_limit\'-parameters are mutually exclusive!'
        )

    # Create permutated variants of 'rj_desc' and 'rj_d_desc'.
    rj_desc_perms = np.reshape(
        np.tile(R_desc[j, :], n_perms)[tril_perms_lin], (n_perms, -1), order='F'
    )

    rj_d_desc = desc_func.d_desc_from_comp(R_d_desc[j, :, :])[0][
        :, keep_idxs_3n
    ]  # convert descriptor back to full representation

    rj_d_desc_perms = np.reshape(
        np.tile(rj_d_desc.T, n_perms)[:, tril_perms_lin], (-1, dim_d, n_perms)
    )

    mat52_base_div = 3 * sig ** 4
    sqrt5 = np.sqrt(5.0)
    sig_pow2 = sig ** 2

    dim_i_keep = rj_d_desc.shape[1]
    diff_ab_outer_perms = np.empty((dim_d, dim_i_keep))
    diff_ab_perms = np.empty((n_perms, dim_d))
    ri_d_desc = np.zeros((1, dim_d, dim_i))  # must be zeros!
    k = np.empty((dim_i, dim_i_keep))

    exploit_sym = True
    #print("exploit_sym = ", exploit_sym)

    for i in range(j,n_train):

        blk_i = slice(i * dim_i, (i + 1) * dim_i)
        #print("blk_i = ", blk_i)

        # diff_ab_perms = R_desc[i, :] - rj_desc_perms
        np.subtract(R_desc[i, :], rj_desc_perms, out=diff_ab_perms)

        norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1)
        mat52_base_perms = np.exp(-norm_ab_perms / sig) / mat52_base_div * 5

        # diff_ab_outer_perms = 5 * np.einsum(
        #    'ki,kj->ij',
        #    diff_ab_perms * mat52_base_perms[:, None],
        #    np.einsum('ik,jki -> ij', diff_ab_perms, rj_d_desc_perms)
        # )
        np.einsum(
            'ki,kj->ij',
            diff_ab_perms * mat52_base_perms[:, None] * 5,
            np.einsum('ki,jik -> kj', diff_ab_perms, rj_d_desc_perms),
            out=diff_ab_outer_perms,
        )

        diff_ab_outer_perms -= np.einsum(
            'ikj,j->ki',
            rj_d_desc_perms,
            (sig_pow2 + sig * norm_ab_perms) * mat52_base_perms,
        )

        # ri_d_desc = desc_func.d_desc_from_comp(R_d_desc[i, :, :])[0]
        desc_func.d_desc_from_comp(R_d_desc[i, :, :], out=ri_d_desc)

        # K[blk_i, blk_j] = ri_d_desc[0].T.dot(diff_ab_outer_perms)
        np.dot(ri_d_desc[0].T, diff_ab_outer_perms, out=k)
        K[blk_i, blk_j] = k

        # print(k - k2)

        if exploit_sym and (
            cols_m_limit is None or i < cols_m_limit
        ):  # this will never be called with 'keep_idxs_3n' set to anything else than [:]
            K[blk_j, blk_i] = K[blk_i, blk_j].T

    if use_E_cstr:
        raise SystemError("Not supported yet")

    return blk_j.stop - blk_j.start


