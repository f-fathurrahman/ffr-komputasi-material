import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "font.size": 11,
    "savefig.dpi": 150
})

np.random.seed(1234)


# Polynomial model + nonlinearity
def true_model(x):
    return (x + 2)*(x + 2)*(x - 3) + 100*np.sin(x) + np.exp(-2*x**2)

# !!! Make sure that this is the actual degree of true_model
TRUE_DEGREE = 3

def generate_data(Ndata=20, noise_amplitude=1.0):
    print(f"Generating {Ndata} data with noise_amplitude = {noise_amplitude}")
    x = np.linspace(-5.0, 5.0, Ndata)
    y = true_model(x) + noise_amplitude*np.random.randn(Ndata)
    return x, y

Ndata = 500
x, y = generate_data(Ndata=Ndata, noise_amplitude=5.0)
X = x[:,np.newaxis]

plt.plot(x, y, label="data", linewidth=0, marker="o")
plt.legend()
plt.savefig("IMG_10_nonlinear_data.png")
plt.show()

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPRegressor

model1 = make_pipeline(
    PolynomialFeatures(TRUE_DEGREE),
    StandardScaler(),
    LinearRegression()
)
cv_score = cross_val_score(model1, X, y, cv=5, scoring="neg_mean_squared_error")
print(f"cv_score = {cv_score}")
print(f"Linear model mean = {cv_score.mean():18.10f} std = {cv_score.std():18.10f}")


model2 = make_pipeline(
    StandardScaler(),
    MLPRegressor(
        hidden_layer_sizes=(64,32), activation="relu",
        solver="adam", alpha=0.0001, batch_size="auto",
        verbose=False, max_iter=10000
    )
)

#model2.fit(X, y)
cv_score = cross_val_score(model2, X, y, cv=5, scoring="neg_mean_squared_error")
print(f"cv_score = {cv_score}")
print(f"MLP model mean = {cv_score.mean():18.10f} std = {cv_score.std():18.10f}")

