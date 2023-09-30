# Linear model, using JAX
# Stochastic gradient descent without batching
#
# Using feature standardization

# Silent TF_CPP messages
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" 
# Not working for JAX ?

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

# Load data
# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")

features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

features = soldata.loc[:, feature_names].values
labels = soldata.Solubility.values

#
# Standardize the features
#
fstd = np.std(features, axis=0)
fmean = np.mean(features, axis=0)
std_features = (features - fmean) / fstd

Ndata = len(labels)
X = np.hstack([
    np.ones((Ndata,1)),
    std_features
])

# fit using numpy least squares method.
# We only need w
w, *_ = np.linalg.lstsq(X, labels, rcond=-1)


#
# Analyzing model performance: using parity plot
#
predicted_labels = np.dot(X,w)
mse_loss = np.mean( (predicted_labels - labels)**2 )
print("mse_loss = ", mse_loss)


plt.clf()
plt.plot([-100, 100], [-100, 100])
plt.scatter(labels, predicted_labels, s=4, alpha=0.7)
plt.xlabel("Measured Solubility $y$")
plt.ylabel("Predicted Solubility $\hat{y}$")
plt.xlim(-13.5, 2)
plt.ylim(-13.5, 2)
plt.gca().set_aspect("equal", "box")
plt.savefig("IMG_main_lstsq_parity_plot.png", dpi=150)
plt.show()

# slice correlation between predict/labels
# from correlation matrix
corr_coef = np.corrcoef(labels, predicted_labels)[0, 1]
# take off-diagonal element
print("corr_coef = ", corr_coef)

res = np.corrcoef(labels, predicted_labels)
print("Full corr_coef = ")
print(res)
