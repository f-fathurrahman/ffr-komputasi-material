# Linear model, using JAX
# Stochastic gradient descent using all data, without batching
# Using feature standardization, see if wee can use larger learning rate

# Silent TF_CPP messages
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" 
# Not working for JAX ?

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp


# Define the model
def linear_model(x, w, b):
    return jnp.dot(x, w) + b

def loss(y, labels):
    return jnp.mean((y - labels)**2)

def loss_wrapper(w, b, data):
    features = data[0]
    labels = data[1]
    y = linear_model(features, w, b)
    return loss(y, labels)

loss_grad = jax.grad(loss_wrapper, (0, 1))

# Load data
# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")

features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

features = soldata.loc[:, feature_names].values
labels = soldata.Solubility.values

feature_dim = features.shape[1]

print("feature_names = ", feature_names)
print("feature_dim = ", feature_dim)

#
# Standardize the features
#
fstd = np.std(features, axis=0)
fmean = np.mean(features, axis=0)
std_features = (features - fmean) / fstd


#
# Initialize model parameters
#
np.random.seed(12345)
w = np.random.normal(size=feature_dim)
b = 0.0

loss_progress = []
eta = 1e-2
batch_size = 32
N = len(labels)  # number of data points
data = (std_features, labels)
# compute how much data fits nicely into a batch
# and drop extra data
new_N = len(labels) // batch_size * batch_size
num_epochs = 100

for epoch in range(num_epochs):
    # use all data (no batching)
    grad = loss_grad(w, b, (std_features, labels))
    w -= eta * grad[0]
    b -= eta * grad[1]
    L = loss_wrapper(w, b, data)
    loss_progress.append(L)
    print("epoch {:5d} loss = {:.5e}".format(epoch, L))

# Print out last loss
print("Last loss = ", loss_wrapper(w, b, data))

# Plot the results
plt.plot(np.arange(len(loss_progress)), loss_progress)
# why multiply with 50? Because we calculate the loss every i modulo 10
plt.xlabel("Step")
plt.yscale("log")
plt.ylabel("Loss")
plt.grid(True)
plt.title("Loss Curve (using std_features)")
plt.savefig("IMG_loss_progress_05.png", dpi=150)
plt.show()


#
# Analyzing model performance: using parity plot
#
predicted_labels = linear_model(std_features, w, b)

plt.clf()
plt.plot([-100, 100], [-100, 100])
plt.scatter(labels, predicted_labels, s=4, alpha=0.7)
plt.xlabel("Measured Solubility $y$")
plt.ylabel("Predicted Solubility $\hat{y}$")
plt.xlim(-13.5, 2)
plt.ylim(-13.5, 2)
plt.gca().set_aspect("equal", "box")
plt.savefig("IMG_main_05_parity_plot.png", dpi=150)
plt.show()

# slice correlation between predict/labels
# from correlation matrix
corr_coef = np.corrcoef(labels, predicted_labels)[0, 1]
# take off-diagonal element
print("corr_coef = ", corr_coef)

res = np.corrcoef(labels, predicted_labels)
print("Full corr_coef = ")
print(res)
