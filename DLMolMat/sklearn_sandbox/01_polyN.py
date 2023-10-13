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
    return (x**3 - 3*x**2 + 0.5*x + 5)/10
    #return (x - 4)*(x - 3)*(x + 2)/20

def generate_data(noise_amplitude=1.0):
    Ndata = 20
    x = np.linspace(-5.0, 5.0, 20)
    y = true_model(x) + noise_amplitude*np.random.randn(Ndata)
    return x, y

x, y = generate_data(noise_amplitude=1.0)

from sklearn.preprocessing import PolynomialFeatures
TRUE_DEGREE = 3
polyN = PolynomialFeatures(TRUE_DEGREE)
X = polyN.fit_transform(x[:,np.newaxis])
print(X)


from sklearn.linear_model import LinearRegression
model = LinearRegression()
model.fit(X, y)


print(model.coef_)
print(model.intercept_)


xdense = np.linspace(-5.0, 5.0, 100)
Xdense = polyN.fit_transform(xdense[:,np.newaxis])
ydense = model.predict(Xdense) # predicted from model

plt.plot(x, y, label="data", linewidth=0, marker="o")
# you also can use scatter
# I use plot to get different color
plt.plot(xdense, ydense, label="model")
plt.legend()
plt.savefig("IMG_01_data.png")
plt.show()


