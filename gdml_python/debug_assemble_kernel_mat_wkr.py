import numpy as np

def debug_assemble_kernel_mat_wkr(
    K, R_desc, R_d_desc, desc_func,
    j, tril_perms_lin, sig
):
    
    use_E_cstr = False
    exploit_sym = True

    n_train, dim_d = R_d_desc.shape[:2]
    n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)
    dim_i = 3 * n_atoms
    n_perms = int(len(tril_perms_lin) / dim_d)

    blk_j = slice(j * dim_i, (j + 1) * dim_i)
    keep_idxs_3n = slice(None)  # same as [:]

    # Create permutated variants of 'rj_desc' and 'rj_d_desc'.
    rj_desc_perms = np.reshape(
        np.tile(R_desc[j, :], n_perms)[tril_perms_lin], (n_perms, -1), order='F'
    )
    # R_desc[j,:] -> rj_desc_perms  (36,) -> (1,36)
    # simply reshaping if n_perms == 1 ?

    rj_d_desc = desc_func.d_desc_from_comp(R_d_desc[j, :, :])[0][:, keep_idxs_3n]
    # convert descriptor back to full representation

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

    for i in range(j, n_train):

        blk_i = slice(i * dim_i, (i + 1) * dim_i)
        #print("blk_i = ", blk_i)

        # diff_ab_perms = R_desc[i, :] - rj_desc_perms
        np.subtract(R_desc[i, :], rj_desc_perms, out=diff_ab_perms)

        norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1)
        mat52_base_perms = np.exp(-norm_ab_perms / sig) / mat52_base_div * 5

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

        if exploit_sym:
            K[blk_j, blk_i] = K[blk_i, blk_j].T

    if use_E_cstr:
        raise SystemError("Not supported yet")

    return

