import numpy as np
np.random.seed(1234)

filename = "DATASET/ethanol_dft.npz"
dataset = np.load(filename)
n_train = 200



T = dataset["E"]
n = n_train
excl_idxs = None

# Freedman-Diaconis rule
h = 2 * np.subtract(*np.percentile(T, [75, 25])) / np.cbrt(n)
print("h = ", h)
n_bins = int(np.ceil((np.max(T) - np.min(T)) / h)) if h > 0 else 1
n_bins = min(
    n_bins, int(n / 2)
)  # Limit number of bins to half of requested subset size.

bins = np.linspace(np.min(T), np.max(T), n_bins, endpoint=False)
idxs = np.digitize(T, bins)

# Exclude restricted indices.
# .... SKIPPED ....

uniq_all, cnts_all = np.unique(idxs, return_counts=True)
print(uniq_all)
print(cnts_all)

# Remove restricted bin.
# .... SKIPPED ....

# Compute reduced bin counts.
reduced_cnts = np.ceil(cnts_all / np.sum(cnts_all, dtype=float) * n).astype(int)
reduced_cnts = np.minimum(
    reduced_cnts, cnts_all
)  # limit reduced_cnts to what is available in cnts_all
print("reduced_cnts = ", reduced_cnts)

# Reduce/increase bin counts to desired total number of points.
reduced_cnts_delta = n - np.sum(reduced_cnts)
print("reduced_cnts_delta = ", reduced_cnts_delta)

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

