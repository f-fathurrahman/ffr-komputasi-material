# From https://jaxopt.github.io/stable/unconstrained.html

import pandas as pd
import numpy as np

import jax
import jax.numpy as jnp
import jaxopt
from jax import random
from jax import grad, jit, vmap

jax.config.update("jax_enable_x64", True)


def load_data():
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

    return std_features, labels # X, y


# A helper function to randomly initialize weights and biases
# for a dense neural network layer
def random_layer_params(m, n, key, scale=1e-2):
    w_key, b_key = random.split(key)
    return scale * random.normal(w_key, (n, m)), scale * random.normal(b_key, (n,))

# Initialize all layers for a fully-connected neural network with sizes "sizes"
def init_network_params(sizes, key):
    keys = random.split(key, len(sizes))
    return [random_layer_params(m, n, k) for m, n, k in zip(sizes[:-1], sizes[1:], keys)]


# apply to one row
def nonlinear_model(params, x):
    activations = x # initialize activation with input
    for w, b in params[:-1]:
        outputs = jnp.dot(w, activations) + b
        activations = jax.nn.sigmoid(outputs)
    # Last output
    final_w, final_b = params[-1]
    outputs = jnp.dot(final_w, activations) + final_b
    return outputs # no activation at last layer (?)


# Make a batched version of the `predict` function
batched_nonlinear_model = vmap(nonlinear_model, in_axes=(None, 0))



# Define objective function
def loss_func(params, X, y):
    residuals = batched_nonlinear_model(params, X).reshape(-1) - y
    return jnp.mean(residuals**2)

# reshape(-1) is important !!!



X, y = load_data()
layer_sizes = [17, 52, 1]
init_params = init_network_params(layer_sizes, random.PRNGKey(0))

solver = jaxopt.LBFGS(fun=loss_func, maxiter=500)
res = solver.run(init_params, X=X, y=y)

# Alternatively, we could have used one of these solvers as well:
# solver = jaxopt.GradientDescent(fun=ridge_reg_objective, maxiter=500)
# solver = jaxopt.ScipyMinimize(fun=ridge_reg_objective, method="L-BFGS-B", maxiter=500)
# solver = jaxopt.NonlinearCG(fun=ridge_reg_objective, method="polak-ribiere", maxiter=500)

import matplotlib.pyplot as plt

y_pred = batched_nonlinear_model(res.params, X)
mse_loss = np.mean( (y - y_pred)**2 )
print("mse_loss = ", mse_loss)

plt.clf()
plt.plot([-100, 100], [-100, 100])
plt.scatter(y, y_pred, s=4, alpha=0.7)
plt.xlabel("Measured Solubility $y$")
plt.ylabel("Predicted Solubility $\hat{y}$")
plt.xlim(-13.5, 2)
plt.ylim(-13.5, 2)
plt.gca().set_aspect("equal", "box")
plt.savefig("IMG_aqsoldb_all_nonlin.png", dpi=150)
plt.show()

