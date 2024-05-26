import numpy as np
from my_draw_strat_sample import my_draw_strat_sample
from my_desc_from_R import my_desc_from_R
from my_desc import Desc

# Prepare input for assembly_kernel_mat
def prepare_data(filename="DATASET/ethanol_dft.npz", n_train=200):

    dataset = np.load(filename)

    #
    # We create task here (manually)
    #
    train_dataset = dataset
    use_sym = False # set to False to use GDML not sGDML
    valid_dataset = dataset
    n_valid = 1000 # not used
    sig = 20 # integer, kernel length scale
    lam = 1e-10

    # Default parameters
    #perms = None
    use_E = True
    use_E_cstr = False

    use_E_cstr = use_E and use_E_cstr
    print("use_E_cstr = ", use_E_cstr)

    n_atoms = train_dataset["R"].shape[1]
    print("n_atoms = ", n_atoms)

    # if "E" in train_dataset:
    print("Will call draw_strat_sample")
    idxs_train = my_draw_strat_sample(train_dataset["E"], n_train)

    # Assuming same md5_train and md5_valid
    excl_idxs = idxs_train

    idxs_valid = my_draw_strat_sample(
        valid_dataset["E"],
        n_valid,
        excl_idxs=excl_idxs,
    )

    R_train = train_dataset["R"][idxs_train, :, :]
    F_train = train_dataset["F"][idxs_train, :, :]
    task = {
        "type": "t",
        "code_version": "0.5.1",
        "dataset_name": train_dataset["name"].astype(str),
        "dataset_theory": train_dataset["theory"].astype(str),
        "z": train_dataset["z"],
        "R_train": R_train,
        "F_train": F_train,
        "idxs_train": idxs_train,
        #"md5_train": md5_train,
        "idxs_valid": idxs_valid,
        #"md5_valid": md5_valid,
        "sig": sig,
        "lam": lam,
        "use_E": use_E,
        "use_E_cstr": use_E_cstr,
        "use_sym": use_sym,
    }

    if use_E:
        print("Using E_train in task")
        task["E_train"] = train_dataset["E"][idxs_train]

    task["perms"] = np.arange(train_dataset["R"].shape[1])[None,:]  # no symmetries

    # from training step

    task = dict(task)  # make mutable (need this?)
    n_train, n_atoms = task["R_train"].shape[:2]

    desc = Desc(n_atoms, max_processes=1)

    n_perms = task["perms"].shape[0]
    tril_perms = np.array([Desc.perm(p) for p in task["perms"]])

    dim_d = desc.dim

    perm_offsets = np.arange(n_perms)[:, None] * dim_d
    tril_perms_lin = (tril_perms + perm_offsets).flatten("F")

    lat_and_inv = None
    R = task['R_train'].reshape(n_train, -1)
    R_desc, R_d_desc = my_desc_from_R(
        desc, R, lat_and_inv=lat_and_inv
    )

    # Generate label vector.
    y = task["F_train"].ravel().copy()
    y_std = np.std(y)
    y /= y_std

    sig = task["sig"]
    lam = task["lam"]
    use_E_cstr = task["use_E_cstr"]

    n_train_, dim_d = R_d_desc.shape[:2]
    assert n_train_ == n_train
    n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)

    return R_desc, R_d_desc, tril_perms_lin, desc, task, y


import pickle
def export_to_pickle(filename="PREPARE_DATA.pkl"):
    ret_vals = prepare_data()
    with open(filename, "wb") as outfile:
        pickle.dump(ret_vals, outfile)
    print(f"Files are saved to {filename}")

def import_from_pickle(filename="PREPARE_DATA.pkl"):
    with open(filename, "rb") as infile:
        in_data = pickle.load(infile)
    return in_data

if __name__ == "__main__":
    export_to_pickle()
