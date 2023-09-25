# Silent TF_CPP messages
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" 
# Not working for JAX ?

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")
plt.rcParams["savefig.dpi"] = 150

np.random.seed(1234)
# Pandas will also use this as random number state

# Load data
# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")

features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]


N = len(soldata)
error = []
error_std = []
for k in range(2, 25):
    splits = list(range(0, N + N // k, N // k))
    k_error = []
    for i in range(k):
        # slice out segments
        test = soldata[splits[i] : splits[i + 1]]
        test_x, test_y = test[feature_names].values, test["Solubility"].values
        train = pd.concat([soldata[splits[i] :], soldata[splits[i + 1] :]])
        x, y = train[feature_names].values, train["Solubility"].values
        # compute coefficients
        w, *_ = np.linalg.lstsq(x, y, rcond=-1)
        # compute intercept (b)
        b = np.mean(y - np.dot(x, w))
        # compute test error
        k_error.append(np.mean((np.dot(test_x, w) + b - test_y) ** 2))
    error.append(np.mean(k_error))
    error_std.append(np.std(k_error, ddof=1))
plt.errorbar(range(2, 25), error, yerr=error_std, capsize=6)
plt.xlabel("k")
plt.ylabel("Test Error")
plt.title("k-fold cross-validation of soldata")
plt.savefig("IMG_kfold_02.png")
plt.show()