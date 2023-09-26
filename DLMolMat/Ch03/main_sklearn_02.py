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

print("feature_nams = ", feature_names)

# convert from pandas dataframe to numpy arrays
X_all = soldata[feature_names].values
y_all = soldata["Solubility"].values

# standardize the features
X_mean = np.mean(X_all, axis=0)
X_std = np.std(X_all, axis=0, ddof=1)
X_all -= X_mean
X_all /= X_std

test_size = int(0.2*X_all.shape[0])
print("test_size = ", test_size)

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression

X, test_X, y, test_y = train_test_split(X_all, y_all, test_size=test_size)
print("train size: ", X.shape[0])

model = LinearRegression()
model.fit(X, y)

R2_score = model.score(X, y)
print("R2_score on training data: ", R2_score)

y_pred = model.predict(X)
RMSE = np.sqrt(np.mean((y_pred - y)**2))
print("RMSE training: ", RMSE)

y_pred = model.predict(test_X)
RMSE = np.sqrt(np.mean((y_pred - test_y)**2))
print("RMSE test: ", RMSE)
