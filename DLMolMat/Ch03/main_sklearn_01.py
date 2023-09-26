import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

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

# Get 50 points and split into train/test
sample = soldata.sample(50, replace=False)
train = sample[:25]
test = sample[25:]


def plot_data(x, y, test_x, test_y, feature_names):
    Nfeatures = len(feature_names)
    for i in range(Nfeatures):
        plt.clf()
        plt.scatter(x[:,i], y[:], label="train")
        plt.scatter(test_x[:,i], test_y[:], label="test")
        plt.xlabel(feature_names[i] + " (standardized)")
        plt.ylabel("Solubility")
        plt.legend()
        plt.grid(True)
        filename = "IMG_sample50_" + feature_names[i] + ".png"
        plt.savefig(filename)

# convert from pandas dataframe to numpy arrays
x = train[feature_names].values
y = train["Solubility"].values
test_x = test[feature_names].values
test_y = test["Solubility"].values


# standardize the features using only train
x_mean = np.mean(x, axis=0)
x_std = np.std(x, axis=0, ddof=1)
#
test_x -= x_mean
test_x /= x_std
#
x -= x_mean
x /= x_std

plot_data(x, y, test_x, test_y, feature_names)


# Should be close to 0 and 1
print("x mean = ", np.mean(x, axis=0))
print("x std  = ", np.std(x, axis=0, ddof=1))

from sklearn.linear_model import LinearRegression
model = LinearRegression()
model.fit(x, y)

R2_score = model.score(x, y)
print("R2_score on training data: ", R2_score)

y_pred = model.predict(x)
RMSE = np.sqrt(np.mean((y_pred - y)**2))
print("RMSE training: ", RMSE)

y_pred = model.predict(test_x)
RMSE = np.sqrt(np.mean((y_pred - test_y)**2))
print("RMSE test: ", RMSE)
