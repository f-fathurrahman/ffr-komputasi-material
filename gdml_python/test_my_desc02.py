import numpy as np

filename = "DATASET/ethanol_dft.npz"
dataset = np.load(filename)

z = dataset["z"] # not changed
idx_data = range(5)
R = dataset["R"][idx_data]
F = dataset["E"][idx_data]
F = dataset["F"][idx_data]

n_atoms = len(z)
desc_dim = n_atoms * (n_atoms - 1)//2

from my_desc_from_R import my_desc_from_R, d_desc_from_comp
lat_and_inv = None
n_train = len(idx_data)
R_lin = R.reshape(n_train, -1) # it is linear?
R_desc, R_d_desc = my_desc_from_R(
    desc_dim, R_lin, lat_and_inv=lat_and_inv
)

# This is special case (?) for n_perms = 1
n_perms = 1
tril_perms_lin = np.arange(0, 36, dtype=np.int32)
j = 1
keep_idxs_3n = slice(None)  # same as [:]
# convert descriptor back to full representation
rj_d_desc = d_desc_from_comp(n_atoms, R_d_desc[j,:, :])[0][:, keep_idxs_3n]
# here R_d_desc is given only for one data, current row/column index j
rj_d_desc_perms = np.reshape(
    np.tile(rj_d_desc.T, n_perms)[:, tril_perms_lin], (-1, desc_dim, n_perms)
)

print("Compressed R_d_desc.shape = ", R_d_desc.shape)
print("Uncompressed rj_d_desc.shape = ", rj_d_desc.shape)
print("Uncompressed rj_d_desc_perms.shape = ", rj_d_desc_perms.shape)

