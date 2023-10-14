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
    #return 0.1*(x + 2)*(x + 2)*(x - 3) + 10*np.sin(x) + 10*np.exp(-2*x**2)
    #return 10*np.exp(-0.5*x**2)
    return np.sin(2*x)/(2*x)

# !!! Make sure that this is the actual degree of true_model
TRUE_DEGREE = 10

def generate_data(Ndata=20, noise_amplitude=1.0):
    print(f"Generating {Ndata} data with noise_amplitude = {noise_amplitude}")
    x = np.linspace(-5.0, 5.0, Ndata)
    y = true_model(x) + noise_amplitude*np.random.randn(Ndata)
    return x, y

Ndata = 500
x, y = generate_data(Ndata=Ndata, noise_amplitude=0.1)
X = x[:,np.newaxis]

plt.plot(x, y, label="data", linewidth=0, marker="o")
plt.legend()
plt.savefig("IMG_10_nonlinear_data.png")
plt.show()

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor

# We don't make any dataset splitting

model1 = make_pipeline(
    PolynomialFeatures(TRUE_DEGREE),
    StandardScaler(),
    LinearRegression()
)
model1.fit(X, y)
R2_score = model1.score(X, y)
print("Linear model R2 = ", R2_score)

ypred = model1.predict(X)
plt.clf()
plt.scatter(y, ypred)
plt.xlabel("Data")
plt.ylabel("Prediction$")
plt.gca().set_aspect("equal", "box")
plt.title("Linear model")
plt.savefig("IMG_11_parity_plot_linear.png", dpi=150)
plt.show()

plt.clf()
plt.scatter(x, ypred, label="prediction", alpha=0.5)
plt.scatter(x, y, label="data", alpha=0.5)
plt.legend()
plt.title("Linear model")
plt.savefig("IMG_11_plot_linear.png", dpi=150)
plt.show()

"""
model2 = make_pipeline(
    StandardScaler(),
    MLPRegressor(
        hidden_layer_sizes=(64,32), activation="relu",
        solver="adam", alpha=0.0001, batch_size="auto",
        verbose=False, max_iter=10000
    )
)
model2.fit(X, y)
R2_score = model2.score(X, y)
print("Nonlinear model R2 = ", R2_score)

ypred = model2.predict(X)
plt.clf()
plt.scatter(y, ypred)
plt.xlabel("Data")
plt.ylabel("Prediction$")
plt.gca().set_aspect("equal", "box")
plt.title("Nonlinear model")
plt.savefig("IMG_11_parity_plot_nonlinear.png", dpi=150)
plt.show()


plt.clf()
plt.scatter(x, ypred, label="prediction", alpha=0.5)
plt.scatter(x, y, label="data", alpha=0.5)
plt.legend()
plt.title("Nonlinear model")
plt.savefig("IMG_11_plot_nonlinear.png", dpi=150)
plt.show()
"""

