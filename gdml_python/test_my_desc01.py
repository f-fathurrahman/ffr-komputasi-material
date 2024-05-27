import numpy as np

filename = "DATASET/ethanol_dft.npz"
dataset = np.load(filename)

z = dataset["z"] # not changed
idx_data = 0
R = dataset["R"][idx_data]
F = dataset["E"][idx_data]
F = dataset["F"][idx_data]

n_atoms = len(z)
desc_dim = n_atoms * (n_atoms - 1)//2

from my_desc_from_R import my_desc_from_R, d_desc_from_comp
lat_and_inv = None
n_train = 1
R_lin = R.reshape(n_train, -1) # it is linear?
R_desc, R_d_desc = my_desc_from_R(
    desc_dim, R_lin, lat_and_inv=lat_and_inv
)
