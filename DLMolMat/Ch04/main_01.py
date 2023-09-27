import pandas as pd
import mordred, mordred.descriptors
import rdkit

import numpy as np
import jax
import jax.numpy as jnp

toxdata = pd.read_csv("../DATASET/clintox.csv")
print(toxdata.head())

calc = mordred.Calculator(mordred.descriptors, ignore_3D=True)
# Samples
molecules = [rdkit.Chem.MolFromSmiles(smi) for smi in toxdata.smiles]

# the invalid molecules were None, so we'll just
# use the fact the None is False in Python
valid_mol_idx = [bool(m) for m in molecules]
valid_mols = [m for m in molecules if m]
print("Number of valid molecules: ", len(valid_mols))

# Requires pandas 1.x
#print("Calculating features:") # This takes quite long time
#features = calc.pandas(valid_mols, nproc=1)
#print("Finished calculating features")

#features = pd.read_csv("../DATASET/features_mordred.csv", low_memory=True)
features = pd.read_pickle("../DATASET/features_mordred.pkl")

features.dropna(inplace=True, axis=1)

labels = toxdata[valid_mol_idx].FDA_APPROVED
features -= features.mean(numeric_only=True)
features /= features.std(numeric_only=True)

# we have some nans in features, likely because std was 0
features.dropna(inplace=True, axis=1)
print(f"We have {len(features.columns)} features per molecule")


def perceptron(x, w, b):
    v = jnp.dot(x, w) + b
    y = jnp.where(v > 0, x=jnp.ones_like(v), y=jnp.zeros_like(v))
    return y



def loss(y, yhat):
    return jnp.mean(jnp.abs(y - yhat))


def loss_wrapper(w, b, x, y):
    yhat = perceptron(x, w, b)
    return loss(y, yhat)


loss_grad = jax.grad(loss_wrapper, (0, 1))


batch_size = 32
train_N = int(len(labels) * 0.8)


N = len(labels)
batch_idx = range(0, train_N, batch_size)
w = np.random.normal(size=len(features.columns))
b = 0.0

loss_grad = jax.grad(loss_wrapper, (0, 1))


test_x = features[train_N:].values.astype(np.float32)
test_y = labels[train_N:].values


loss_grad(w, b, test_x, test_y)


def bin_classifier(x, w, b):
    v = jnp.dot(x, w) + b
    y = jax.nn.sigmoid(v)
    return y


def cross_ent(y, yhat):
    return jnp.mean(-(y * jnp.log(yhat + 1e-10) + (1 - y) * jnp.log(1 - yhat + 1e-10)))


def loss_wrapper(w, b, x, y):
    yhat = bin_classifier(x, w, b)
    return cross_ent(y, yhat)


loss_grad = jax.grad(loss_wrapper, (0, 1))
w = np.random.normal(scale=0.01, size=len(features.columns))
b = 1.0

loss_progress = []
test_loss_progress = []
eta = 0.2
for epoch in range(5):
    for i in range(len(batch_idx) - 1):
        x = features[batch_idx[i] : batch_idx[i + 1]].values.astype(np.float32)
        y = labels[batch_idx[i] : batch_idx[i + 1]].values
        grad = loss_grad(w, b, x, y)
        w -= eta * grad[0]
        b -= eta * grad[1]
        loss_progress.append(loss_wrapper(w, b, x, y))
        test_loss_progress.append(loss_wrapper(w, b, test_x, test_y))

import matplotlib.pyplot as plt

plt.clf()
plt.plot(loss_progress, label="Training Loss")
plt.plot(test_loss_progress, label="Testing Loss")
plt.xlabel("Step")
plt.legend()
plt.ylabel("Loss")
plt.show()

