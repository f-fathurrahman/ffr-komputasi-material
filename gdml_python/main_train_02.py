import numpy as np

import numpy as np
from my_sgdml.utils.desc import Desc
from my_sgdml.train import GDMLTrain

np.random.seed(1234)

dataset = np.load("DATASET/ethanol_dft.npz")
#dataset = np.load("DATASET/benzene2017_dft.npz")
for k in dataset.keys():
    print("key = ", k)
n_train = 200

gdml_train = GDMLTrain(max_processes=1)
# global variable glob is here?

print("Create task")
task = gdml_train.create_task(dataset, n_train,
        valid_dataset=dataset, n_valid=1000,
        sig=20, lam=1e-10, use_sym=False)
print("End create task")


n_train, n_atoms = task['R_train'].shape[:2]
desc = Desc(n_atoms, max_processes=1)

n_perms = task['perms'].shape[0]
tril_perms = np.array([Desc.perm(p) for p in task['perms']])
# Array of size N*(N-1)/2 containing the corresponding descriptor permutation.
# symmetric matrix, lower triangular array

# Convert to flat array
dim_d = desc.dim
perm_offsets = np.arange(n_perms)[:, None] * dim_d
tril_perms_lin = (tril_perms + perm_offsets).flatten('F')

#
# Descriptor is calculated here
#

from my_desc_from_R import my_desc_from_R
lat_and_inv = None
R = task['R_train'].reshape(n_train, -1)
desc_dim = (n_atoms * (n_atoms - 1)) // 2
R_desc, R_d_desc = my_desc_from_R(
    desc_dim, R, lat_and_inv=lat_and_inv
)

# Generate label vector.
E_train_mean = None
y = task['F_train'].ravel().copy() # the original F_train is not modified
# Need this?
if task['use_E'] and task['use_E_cstr']:
    print("use_E is True")
    E_train = task['E_train'].ravel().copy()
    E_train_mean = np.mean(E_train)
    y = np.hstack((y, -E_train + E_train_mean))
y_std = np.std(y)
print("y_std = ", y_std)
y /= y_std


print("Using analytic solver")
#analytic = Analytic(gdml_train, desc)
#alphas = analytic.solve(task, R_desc, R_d_desc, tril_perms_lin, y)
from my_analytic_solve import my_analytic_solve
alphas = my_analytic_solve(desc, task, R_desc, R_d_desc, tril_perms_lin, y)
print("End of finding parameters")

print(type(alphas))
print("alphas.shape = ", alphas.shape)
print("average alpha = ", np.average(alphas))

alphas_E = None
alphas_F = alphas
if task['use_E_cstr']:
    alphas_E = alphas[-n_train:]
    alphas_F = alphas[:-n_train]


solver_keys = {}
print("After finding the parameters, create model which will be returned")
print("solver_keys = ", solver_keys)

model = gdml_train.create_model(
    task,
    'analytic',
    R_desc,
    R_d_desc,
    tril_perms_lin,
    y_std,
    alphas_F,
    alphas_E=alphas_E,
)
model.update(solver_keys) # not really used
# model is a dict

# Recover integration constant.
# Note: if energy constraints are included in the kernel (via 'use_E_cstr'), do not
# compute the integration constant, but simply set it to the mean of the training energies
# (which was subtracted from the labels before training).
if model['use_E']:
    c = (
        gdml_train._recov_int_const(model, task, R_desc=R_desc, R_d_desc=R_d_desc)
        if E_train_mean is None
        else E_train_mean
    )
    print("Recover integration constant: c = ", c)
    if c is None:
        # Something does not seem right. Turn off energy predictions for this model, only output force predictions.
        model['use_E'] = False
    else:
        model['c'] = c


