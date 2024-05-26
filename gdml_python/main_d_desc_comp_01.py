import numpy as np
import numpy.random
numpy.random.seed(1234)

import matplotlib.style
matplotlib.style.use("dark_background")


def my_d_desc_from_comp(desc, R_d_desc, out=None):

    if R_d_desc.ndim == 2:
        R_d_desc = R_d_desc[None, ...]

    n = R_d_desc.shape[0]
    i, j = desc.tril_indices
    # desc.tril_indices are lower triangular indices, i and j indices for atoms?

    if out is None:
        out = np.zeros((n, desc.dim, desc.n_atoms, 3))
    else:
        out = out.reshape(n, desc.dim, desc.n_atoms, 3)

    out[:, desc.dim_range, j, :] = R_d_desc
    out[:, desc.dim_range, i, :] = -R_d_desc

    return out.reshape(-1, desc.dim, desc.dim_i)

# Load data
from prepare_data import import_from_pickle
R_desc, R_d_desc, tril_perms_lin, desc, task, y = import_from_pickle()

# Some dimensions or array size
n_train, dim_d = R_d_desc.shape[:2]
n_perms = int(len(tril_perms_lin) / dim_d)

j = 1 # row index
rj_desc_perms = np.reshape(
    np.tile(R_desc[j, :], n_perms)[tril_perms_lin], (n_perms, -1), order='F'
)
# This will simply convert to row matrix (for n_perms==1)

keep_idxs_3n = slice(None)  # same as [:]
# convert R_d_desc to full form
rj_d_desc = my_d_desc_from_comp(desc, R_d_desc[j, :, :])[0][:, keep_idxs_3n]
rj_d_desc_perms = np.reshape(
    np.tile(rj_d_desc.T, n_perms)[:, tril_perms_lin], (-1, dim_d, n_perms)
)
# How `R_d_desc[j,:,:]` become `rj_d_desc` ?


# desc.dim_i = 3 * n_atoms
# # Size of the resulting descriptor vector.
# desc.dim = (n_atoms * (n_atoms - 1)) // 2

# Some results:
#
# In  : R_desc[j,:].shape
# Out : (36,)
#
# In  : rj_desc_perms.shape
# Out : (1, 36)  <-- transposed, add singleton dimension (first dim)
#
# In  : rj_d_desc.shape
# Out : (36, 27)
#
# In  : rj_d_desc_perms.shape
# Out : (27, 36, 1)   <-- transposed, add singleton dimension (last dim)
# 
