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

# Initialize parameters
np.random.seed(12345)
w = np.random.normal(size=feature_dim)
b = 0.0

y = linear_model(features, w, b)
print("Current loss = ", loss(y, labels))


# compute gradients
def loss_wrapper(w, b, data):
    features = data[0]
    labels = data[1]
    y = linear_model(features, w, b)
    return loss(y, labels)

loss_grad = jax.grad(loss_wrapper, (0, 1))

# test it out
print("loss grad = ", loss_grad(w, b, (features, labels)))


loss_progress = []
eta = 1e-6
data = (features, labels)
for i in range(10):
    grad = loss_grad(w, b, data)
    w -= eta * grad[0]
    b -= eta * grad[1]
    loss_progress.append(loss_wrapper(w, b, data))

print("Last loss = ", loss_wrapper(w, b, data))

plt.plot(loss_progress, marker="o")
plt.xlabel("Step")
plt.yscale("log")
plt.ylabel("Loss")
plt.grid(True)
plt.title("Full Dataset Training Curve")
plt.savefig("IMG_loss_progress_01.png", dpi=150)
plt.show()
