import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "font.size": 11,
    "savefig.dpi": 150
})

np.random.seed(1234)

GBL_a1 = 1.2
GBL_a2 = 9.9
GBL_b = 3.1

def true_model(x1, x2):
    return GBL_a1*x1 + GBL_a2*x2 + GBL_b

def generate_data(noise_amplitude=1.0):
    print("Generating data with noise_amplitude = ", noise_amplitude)
    Ndata = 10
    x1 = np.linspace(-5.0, 5.0, Ndata)
    x2 = np.linspace(-5.0, 5.0, Ndata) + np.random.randn(Ndata)
    y = true_model(x1, x2) + noise_amplitude*np.random.randn(Ndata)
    return x1, x2, y




x1, x2, y = generate_data(noise_amplitude=0.2)

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression

X = np.vstack([x1,x2]).T
GG = X.T @ X
print("Condition number of Gramian matrix: ", np.linalg.cond(GG))

model = make_pipeline(
    StandardScaler(),
    LinearRegression()
)
model.fit(X, y)

ypred = model.predict(X)
mse = np.mean((y - ypred)**2)
print("mse = ", mse)
print("R2 = ", model.score(X, y))


