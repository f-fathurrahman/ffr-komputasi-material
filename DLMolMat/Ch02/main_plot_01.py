import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

def plot_data(x, y, feature_names):
    Nfeatures = len(feature_names)
    for i in range(Nfeatures):
        plt.clf()
        plt.scatter(x[:,i], y[:])
        plt.xlabel(feature_names[i] + "(original)")
        plt.ylabel("Solubility")
        plt.grid(True)
        filename = "IMG_soldata_" + feature_names[i] + ".png"
        plt.savefig(filename)

plt.rcParams.update({
    "font.size": 12,
    "savefig.dpi": 150
})

np.random.seed(1234)
# Pandas will also use this as random number state

# Load data
# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")

features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

# convert from pandas dataframe to numpy arrays
x = soldata[feature_names].values
y = soldata["Solubility"].values

plot_data(x, y, feature_names)

