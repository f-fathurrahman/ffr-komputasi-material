# Silent TF_CPP messages
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3" 
# Not working for JAX ?

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp
from jax.example_libraries import optimizers


# Load data
# Original: https://dataverse.harvard.edu/api/access/datafile/3407241
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")

features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

# Get 50 points and split into train/test
sample = soldata.sample(50, replace=False)
train = sample[:25]
test = sample[25:]

# standardize the features using only train
test[feature_names] -= train[feature_names].mean()
test[feature_names] /= train[feature_names].std()
train[feature_names] -= train[feature_names].mean()
train[feature_names] /= train[feature_names].std()

# convert from pandas dataframe to numpy arrays
x = train[feature_names].values
y = train["Solubility"].values
test_x = test[feature_names].values
test_y = test["Solubility"].values


# define our loss function
@jax.jit
def loss(w, b, x, y):
    return jnp.mean((y - jnp.dot(x, w) - b) ** 2)


loss_grad = jax.grad(loss, (0, 1))
w = np.random.normal(size=x.shape[1])
b = 0.0
loss_grad(w, b, x, y)

loss_progress = []
test_loss_progress = []
eta = 0.05
for i in range(2000):
    grad = loss_grad(w, b, x, y)
    w -= eta * grad[0]
    b -= eta * grad[1]
    loss_progress.append(loss(w, b, x, y))
    test_loss_progress.append(loss(w, b, test_x, test_y))
plt.plot(loss_progress, label="Training Loss")
plt.plot(test_loss_progress, label="Testing Loss")

plt.xlabel("Step")
plt.yscale("log")
plt.legend()
plt.ylabel("Loss")
plt.savefig("IMG_01.png", dpi=150)

