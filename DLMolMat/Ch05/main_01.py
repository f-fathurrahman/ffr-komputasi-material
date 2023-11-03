import pandas as pd
import numpy as np
import jax.numpy as jnp
import jax

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "font.size": 12,
    "savefig.dpi": 150
})

# For reproducibility
np.random.seed(0)


# Read data
soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

# Split into train(80) and test(20)
N = len(soldata)
split = int(N * 0.8)
shuffled = soldata.sample(N, replace=False)
train = shuffled[:split]
test = shuffled[split:]

# standardize the features using only train
test[feature_names] -= train[feature_names].mean()
test[feature_names] /= train[feature_names].std()
train[feature_names] -= train[feature_names].mean()
train[feature_names] /= train[feature_names].std()


# Kernel definition
def kernel(x1, x2):
    return jnp.sqrt(jnp.mean((x1 - x2) ** 2))


# Model definition
def model(x, train_x, w, b):
    # make vectorized version of kernel
    vkernel = jax.vmap(kernel, in_axes=(None, 0), out_axes=0)
    # compute kernel with all training data
    s = vkernel(x, train_x)
    # dual form
    yhat = jnp.dot(s, w) + b
    return yhat


# make batched version that can handle multiple xs
batch_model = jax.vmap(model, in_axes=(0, None, None, None), out_axes=0)


# The loss function
@jax.jit
def loss(w, b, train_x, x, y):
    return jnp.mean((batch_model(x, train_x, w, b) - y) ** 2)

# Find the derivative
loss_grad = jax.grad(loss, (0, 1))

# convert from pandas dataframe to numpy arrays
train_x = train[feature_names].values
train_y = train["Solubility"].values
test_x = test[feature_names].values
test_y = test["Solubility"].values

# Parameter for training
eta = 1e-5
batch_size = 32
epochs = 10

# reshape into batches
batch_num = train_x.shape[0] // batch_size
# first truncate data so it's whole nubmer of batches
trunc = batch_num * batch_size
train_x = train_x[:trunc]
train_y = train_y[:trunc]
# split into batches
x_batches = train_x.reshape(-1, batch_size, train_x.shape[-1])
y_batches = train_y.reshape(-1, batch_size)


# make trainable parameters
# w = np.random.normal(scale = 1e-30, size=train_x.shape[0])
w = np.zeros(train_x.shape[0])
b = np.mean(train_y)  # just set to mean, since it's a good first guess



loss_progress = []
test_loss_progress = []
for _ in range(epochs):
    # go in random order
    for i in np.random.randint(0, batch_num - 1, size=batch_num):
        # update step
        x = x_batches[i]
        y = y_batches[i]
        loss_progress.append(loss(w, b, train_x, x, y))
        test_loss_progress.append(loss(w, b, train_x, test_x, test_y))
        grad = loss_grad(w, b, train_x, x, y)
        w -= eta * grad[0]
        b -= eta * grad[1]
plt.plot(loss_progress, label="Training Loss")
plt.plot(test_loss_progress, label="Testing Loss")

plt.xlabel("Step")
plt.yscale("log")
plt.legend()
plt.ylabel("Loss")
plt.show()


