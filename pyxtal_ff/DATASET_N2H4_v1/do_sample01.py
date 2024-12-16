import pickle
import numpy as np

# Adapted from sGDML
def draw_strat_sample(T, n, excl_idxs=None):

    if excl_idxs is None or len(excl_idxs) == 0:
        excl_idxs = None

    if n == 0:
        return np.array([], dtype=np.uint)

    if T.size == n:  # TODO: this only works if excl_idxs=None
        assert excl_idxs is None
        return np.arange(n)

    if n == 1:
        idxs_all_non_excl = np.setdiff1d(
            np.arange(T.size), excl_idxs, assume_unique=True
        )
        return np.array([np.random.choice(idxs_all_non_excl)])

    # Freedman-Diaconis rule
    h = 2 * np.subtract(*np.percentile(T, [75, 25])) / np.cbrt(n)
    n_bins = int(np.ceil((np.max(T) - np.min(T)) / h)) if h > 0 else 1
    n_bins = min(
        n_bins, int(n / 2)
    )  # Limit number of bins to half of requested subset size.

    bins = np.linspace(np.min(T), np.max(T), n_bins, endpoint=False)
    idxs = np.digitize(T, bins)

    # Exclude restricted indices.
    if excl_idxs is not None and excl_idxs.size > 0:
        idxs[excl_idxs] = n_bins + 1  # Impossible bin.

    uniq_all, cnts_all = np.unique(idxs, return_counts=True)

    # Remove restricted bin.
    if excl_idxs is not None and excl_idxs.size > 0:
        excl_bin_idx = np.where(uniq_all == n_bins + 1)
        cnts_all = np.delete(cnts_all, excl_bin_idx)
        uniq_all = np.delete(uniq_all, excl_bin_idx)

    # Compute reduced bin counts.
    reduced_cnts = np.ceil(cnts_all / np.sum(cnts_all, dtype=float) * n).astype(int)
    reduced_cnts = np.minimum(
        reduced_cnts, cnts_all
    )  # limit reduced_cnts to what is available in cnts_all

    # Reduce/increase bin counts to desired total number of points.
    reduced_cnts_delta = n - np.sum(reduced_cnts)

    while np.abs(reduced_cnts_delta) > 0:

        # How many members can we remove from an arbitrary bucket
        # without any bucket with more than one member going to zero?
        max_bin_reduction = np.min(reduced_cnts[np.where(reduced_cnts > 1)]) - 1

        # Generate additional bin members to fill up/drain bucket counts
        # of subset. This array contains (repeated) bucket IDs.
        outstanding = np.random.choice(
            uniq_all,
            min(max_bin_reduction, np.abs(reduced_cnts_delta)),
            p=(reduced_cnts - 1) / np.sum(reduced_cnts - 1, dtype=float),
            replace=True,
        )
        uniq_outstanding, cnts_outstanding = np.unique(
            outstanding, return_counts=True
        )  # Aggregate bucket IDs.

        outstanding_bucket_idx = np.where(
            np.in1d(uniq_all, uniq_outstanding, assume_unique=True)
        )[
            0
        ]  # Bucket IDs to Idxs.
        reduced_cnts[outstanding_bucket_idx] += (
            np.sign(reduced_cnts_delta) * cnts_outstanding
        )
        reduced_cnts_delta = n - np.sum(reduced_cnts)

    # Draw examples for each bin.
    idxs_train = np.empty((0,), dtype=int)
    for uniq_idx, bin_cnt in zip(uniq_all, reduced_cnts):
        idx_in_bin_all = np.where(idxs.ravel() == uniq_idx)[0]
        idxs_train = np.append(
            idxs_train, np.random.choice(idx_in_bin_all, bin_cnt, replace=False)
        )

    return idxs_train




np.random.seed(1234)

with open("ALL_ATOMS.pkl", "rb") as f:
    all_atoms = pickle.load(f)

Ndata = len(all_atoms)
Natoms = len(all_atoms[0])

energies = np.zeros((Ndata))
for i in range(Ndata):
    energies[i] = all_atoms[i].calc.get_potential_energy()

# Number of training data
Ntrain = 200

# Get idx train
idx_train = draw_strat_sample(energies, Ntrain)
energies_train = energies[idx_train]

# Prepare for training data
atoms_train = []
for i in idx_train:
    atoms_train.append(all_atoms[i])

# write as pickle file
#with open("TEMP_ATOMS_TRAIN.pkl", "wb") as f:
#    pickle.dump(atoms_train, f)


filename_xyz = "TEMP_ATOMS_TRAIN.xyz"
atoms_train[0].write(filename_xyz)
for i in range(1,len(atoms_train)):
    atoms_train[i].write(filename_xyz, append=True)

