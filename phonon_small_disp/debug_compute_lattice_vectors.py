import numpy as np
#supercell = [9, 1, 1]
supercell = [3, 4, 2]
offset = 0

# Lattice vectors relevative to the reference cell
R_cN = np.indices(supercell).reshape(3, -1)
N_c = np.array(supercell)[:, np.newaxis]
if offset == 0:
    R_cN += N_c // 2
    R_cN %= N_c
R_cN -= N_c // 2
