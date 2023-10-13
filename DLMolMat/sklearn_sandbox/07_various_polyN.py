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

def generate_data(noise_amplitude=1.0):
    print("Generating data with noise_amplitude = ", noise_amplitude)
    Ndata = 20
    x = np.linspace(-5.0, 5.0, 20)
    y = true_model(x) + noise_amplitude*np.random.randn(Ndata)
    return x, y

x, y = generate_data(noise_amplitude=5.0)

plt.plot(x, y, label="data", linewidth=0, marker="o")
plt.legend()
plt.savefig("IMG_01_data.png")
plt.show()


from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression


# Plot data first
plt.plot(x, y, label="data", linewidth=0, marker="o")

for deg in range(1,10):
    model = make_pipeline(
        PolynomialFeatures(deg),
        StandardScaler(),
        LinearRegression()
    )
    model.fit(x[:,np.newaxis], y)

    ypred = model.predict(x[:,np.newaxis])
    mse = np.mean((y - ypred)**2)
    R2 = model.score(x[:,np.newaxis], y)
    print(f"deg = {deg:3d} mse = {mse:18.10e} R2 = {R2:18.10f}")

    xdense = np.linspace(-6.0, 6.0, 100)
    ydense = model.predict(xdense[:,np.newaxis]) # predicted from model
    plt.plot(xdense, ydense, label=f"degree-{deg}")

plt.legend()
plt.title("06 pipeline")
plt.savefig("IMG_07_data.png")
plt.show()


