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
    return -x**3 + x**2 + x + 2
    #return -x**10 + 2
    #return x**20 - 10*x**9 + 1

# !!! Make sure that this is the actual degree of true_model
TRUE_DEGREE = 3

def generate_data(noise_amplitude=1.0):
    print("Generating data with noise_amplitude = ", noise_amplitude)
    Ndata = 20
    x = np.linspace(-5.0, 5.0, 20)
    y = true_model(x) + noise_amplitude*np.random.randn(Ndata)
    return x, y

x, y = generate_data(noise_amplitude=10.0)

plt.plot(x, y, label="data", linewidth=0, marker="o")
plt.legend()
plt.savefig("IMG_01_data.png")
plt.show()


from sklearn.preprocessing import PolynomialFeatures, StandardScaler

polyN = PolynomialFeatures(TRUE_DEGREE)
polyN.fit(x[:,np.newaxis])
X = polyN.transform(x[:,np.newaxis])

std_scaler = StandardScaler()
std_scaler.fit(X)
Xstd = std_scaler.transform(X)

from sklearn.linear_model import LinearRegression
model = LinearRegression()
model.fit(Xstd, y)

# These will be different from true model parameters because
# we use feature transformation
print("Coefficients = ")
print(model.coef_)
print("Intercept = ")
print(model.intercept_)

ypred = model.predict(Xstd)
mse = np.mean((y - ypred)**2)
print("mse = ", mse)
print("R2 = ", model.score(Xstd, y))

xdense = np.linspace(-5.0, 5.0, 100)

# Do the same preprocessing or feature transformations
# We do not create new preprocessors
Xdense = polyN.transform(xdense[:,np.newaxis])
Xdense = std_scaler.transform(Xdense)

# This is wrong !!!
#Xdense = PolynomialFeatures(TRUE_DEGREE).fit_transform(xdense[:,np.newaxis])
#Xdense = StandardScaler().fit_transform(Xdense)

ydense = model.predict(Xdense) # predicted from model

plt.plot(x, y, label="data", linewidth=0, marker="o")
# you also can use scatter
# I use plot to get different color
plt.plot(xdense, ydense, label="model")
plt.legend()
plt.savefig("IMG_01_data.png")
plt.show()


