import numpy as np
import matplotlib.pyplot as plt

from catlearn.preprocess.scaling import standardize, target_standardize
from catlearn.regression import GaussianProcess
from catlearn.regression.cost_function import get_error

def my_func(x):
    y = x - 50.0
    p = (y + 4) * (y + 4) * (y + 1) * (y - 1) * (y - 3.5) * (y - 2) * (y - 1)
    p = p + 40.0 * y + 80.0 * np.sin(10.0*x)
    return 1.0/20.0 * p + 500

Ntrain = 33
noise_magnitude = 1.0

# Randomly generate training datapoints
train = 7.6 * np.random.sample((Ntrain,1)) - 4.2 + 50

target = np.array(my_func(train))

# Random noise
target = target + noise_magnitude * np.random.randn(Ntrain, 1)

# Generate test datapoints x
Ntest = 513
test = np.vstack( np.linspace(np.min(train) - 0.1, np.max(train) + 0.1, Ntest) )
print(test.shape)

# standard deviations of training data and targets
stdx = np.std(train)
stdy = np.std(target)
tstd = 2.0

# Standardize training and test data
std = standardize(train_matrix=train, test_matrix=test)

# std is a dict with keys:
# dict_keys(['mean', 'std', 'train', 'test'])

# Standardize target data
train_targets = target_standardize(target)

# Predictions will now be made on the standardized scale.

# Store the known underlying function for plotting
linex = np.linspace(np.min(test), np.max(test), Ntest)
liney = my_func(linex)

#
# Example 1: Biased model
#

sdt2 = 0.001
w2 = 0.05 # too small width lead to overfit
kdict = [ {"type": "gaussian", "width": w2} ]

gp = GaussianProcess(kernel_list=kdict, regularization=sdt2,
                     train_fp=std["train"],
                     train_target=train_targets["target"],
                     optimize_hyperparameters=False)

over_fit = gp.predict(test_fp=std["test"], uncertainty=True)

# Scale predictions back to the original scale
over_prediction = np.vstack(over_fit["prediction"]) * train_targets["std"] + train_targets["mean"]
over_certainty = np.vstack(over_fit["uncertainty_with_reg"]) * train_targets["std"]

# Get average errors
error = get_error(over_prediction.reshape(-1), my_func(test).reshape(-1))

# Get confidence interval on predictions
upper = over_prediction + over_certainty*tstd
lower = over_prediction - over_certainty*tstd

plt.figure(0)
plt.plot(linex, liney, "-", lw=1, color="black")
plt.plot(train, target, "o", alpha=0.5, color="black")
plt.plot(test, over_prediction, "b-", lw=1, alpha=0.4)
plt.fill_between(np.hstack(test), np.hstack(upper), np.hstack(lower), 
                 interpolate=True, color="red", alpha=0.2)
plt.title("Biased kernel regression model. \n" +
          "w: {0:.3f}, r: {1:.3f}".format(w2*stdx, sdt2*stdy))
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("tight")

plt.savefig("IMG_overfit.pdf")

