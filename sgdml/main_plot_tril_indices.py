import matplotlib.pyplot as plt
import numpy as np
np.random.seed(1234)

from prepare_data import import_from_pickle
R_desc, R_d_desc, tril_perms_lin, desc, task, y = import_from_pickle()

idx_row, idx_col = desc.tril_indices
A = np.zeros((9,9))
A[idx_row,idx_col] = 1.0
A[idx_col,idx_row] = -1.0
plt.matshow(A)

