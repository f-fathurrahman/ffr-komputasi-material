# From https://jaxopt.github.io/stable/unconstrained.html

import pandas as pd
import numpy as np

import jax
import jax.numpy as jnp
import jaxopt
jax.config.update("jax_enable_x64", True)


def linear_model(X, params):
    W, b = params
    return jnp.dot(X, W) + b


# Define objective function
def ridge_reg_objective(params, l2reg, X, y):
    residuals = linear_model(X, params) - y
    return jnp.mean(residuals ** 2) + 0.5 * l2reg * jnp.sum(params[0]**2)



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



X, y = load_data()
l2reg = 0.0
w_init = np.random.rand(X.shape[1])
b = np.array([0.0])
init_params = (w_init, b)

solver = jaxopt.LBFGS(fun=ridge_reg_objective, maxiter=500)
res = solver.run(init_params, l2reg=l2reg, X=X, y=y)

# Alternatively, we could have used one of these solvers as well:
# solver = jaxopt.GradientDescent(fun=ridge_reg_objective, maxiter=500)
# solver = jaxopt.ScipyMinimize(fun=ridge_reg_objective, method="L-BFGS-B", maxiter=500)
# solver = jaxopt.NonlinearCG(fun=ridge_reg_objective, method="polak-ribiere", maxiter=500)

import matplotlib.pyplot as plt

y_pred = linear_model(X, res.params)
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
plt.savefig("IMG_aqsoldb_all.png", dpi=150)
plt.show()

