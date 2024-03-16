import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

# input: x
# parameter model: w, b
def linear_model(x, w, b):
    return jnp.dot(x, w) + b

# define loss
def loss(y, labels):
    return jnp.mean((y - labels)**2)

# compute gradients
def loss_wrapper(w, b, data):
    features = data[0]
    labels = data[1]
    y = linear_model(features, w, b)
    return loss(y, labels)

loss_grad = jax.grad(loss_wrapper, (0, 1))


soldata = pd.read_csv("../DATASET/curated-solubility-dataset.csv")
features_start_at = list(soldata.columns).index("MolWt")
feature_names = soldata.columns[features_start_at:]

# convert data into features, labels
features = soldata.loc[:, feature_names].values
labels = soldata.Solubility.values

feature_dim = features.shape[1]

# initialize our paramaters
np.random.seed(1234)
w = np.random.normal(size=feature_dim)
b = 0.0

fstd = np.std(features, axis=0) # stdev
fmean = np.mean(features, axis=0) # mean
std_features = (features - fmean) / fstd

loss_progress = []
eta = 1e-2
data = (std_features, labels)
print("INFO: no batch, with feature standardization, learning rate = ", eta)
for i in range(10):
    dLdw, dLdb = loss_grad(w, b, data)
    w -= eta * dLdw
    b -= eta * dLdb
    loss_val = loss_wrapper(w, b, data) # Hitung loss
    print("%5d %18.10f" % (i, loss_val))
    loss_progress.append(loss_val)
plt.plot(loss_progress, marker="o")
plt.savefig("IMG_training_04.pdf", dpi=150)
