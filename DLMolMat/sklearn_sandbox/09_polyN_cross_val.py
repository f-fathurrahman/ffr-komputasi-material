import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "font.size": 11,
    "savefig.dpi": 150
})

np.random.seed(1234)


# Polynomial model
def true_model(x):
    return (x + 2)*(x + 2)*(x - 3)

# !!! Make sure that this is the actual degree of true_model
TRUE_DEGREE = 3

def generate_data(Ndata=20, noise_amplitude=1.0):
    print(f"Generating {Ndata} data with noise_amplitude = {noise_amplitude}")
    x = np.linspace(-5.0, 5.0, Ndata)
    y = true_model(x) + noise_amplitude*np.random.randn(Ndata)
    return x, y

Ndata = 100
x, y = generate_data(Ndata=Ndata, noise_amplitude=1.0)
X = x[:,np.newaxis]

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score

for deg in range(1,10):
    model = make_pipeline(
        PolynomialFeatures(deg),
        StandardScaler(),
        LinearRegression()
    )

    cv_score = cross_val_score(model, X, y, cv=5, scoring="neg_mean_squared_error")
    #cv_score = cross_val_score(model, X, y, cv=5, scoring="r2")
    print(f"deg = {deg:3d} mean = {cv_score.mean():18.10f} std = {cv_score.std():18.10f}")
