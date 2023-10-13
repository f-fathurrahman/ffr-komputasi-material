import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")

plt.rcParams.update({
    "font.size": 11,
    "savefig.dpi": 150
})

#np.random.seed(1234)


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

Ndata = 20
x, y = generate_data(Ndata=Ndata, noise_amplitude=4.0)
X = x[:,np.newaxis]

plt.plot(x, y, label="data", linewidth=0, marker="o")
plt.legend()
plt.savefig("IMG_01_data.png")
plt.show()


from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

# Plot data first
plt.plot(x, y, label="data", linewidth=0, marker="o")


# Split the data
Ndata_test = int(0.2*Ndata)
print("Ndata_test = ", Ndata_test)
Xtrain, Xtest, ytrain, ytest = train_test_split(X, y, test_size=Ndata_test)

losses = {}
for deg in range(1,10):
    model = make_pipeline(
        PolynomialFeatures(deg),
        StandardScaler(),
        LinearRegression()
    )
    model.fit(Xtrain, ytrain)

    # Evaluate performance on test dataset
    ypred = model.predict(Xtest)
    mse = np.mean((ytest - ypred)**2)
    R2 = model.score(Xtest, ytest)
    print(f"deg = {deg:3d} mse = {mse:10.5e} R2 = {R2:10.5f}")

    losses.update({deg : mse})

    # For plotting
    xdense = np.linspace(-5.0, 5.0, 100)
    ydense = model.predict(xdense[:,np.newaxis]) # predicted from model
    plt.plot(xdense, ydense, label=f"degree-{deg}")

# sort the dictionary
losses = sorted(losses.items(), key=lambda x: x[1])

print()
print("Three lowest losses:")
for i in range(3):
    loss = losses[i]
    print(f"{i+1} deg = {loss[0]} mse = {loss[1]}")

plt.legend()
plt.title("06 pipeline")
plt.savefig("IMG_07_data.png")
plt.show()


