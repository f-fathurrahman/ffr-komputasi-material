import numpy as np
import matplotlib.pyplot as plt

from catlearn.preprocess.scaling import standardize, target_standardize
from catlearn.regression import GaussianProcess
from catlearn.regression.cost_function import get_error

def my_func(x):
    A  = 9.1
    ω1 = 0.1
    ω2 = 0.2
    return A + np.sin(ω1*x) + np.cos(ω2*x)

Ntrain = 40
A_noise = 1.0

# Randomly generate training datapoints
#x_train = 7.6 * np.random.sample((Ntrain,1)) - 4.2 + 50
x_train = 5.0*np.random.sample((Ntrain,1)) + 10.0
y_train = my_func(x_train)

# Add noise (normal distributed)
y_train = y_train + A_noise * np.random.randn(Ntrain, 1)

# Generate test datapoints (around x_train)
Ntest = 50
x_test = np.vstack( np.linspace(np.min(x_train)-0.1, np.max(x_train)+0.1, Ntest) )

# standard deviations of training data and targets
std_dev_x = np.std(x_train)
std_dev_y = np.std(y_train)
tstd = 2.0

# Standardize training and test data
std_data = standardize(train_matrix=x_train, test_matrix=x_test)
# In [8]: std_data.keys()                         
# Out[8]: dict_keys(['mean', 'std', 'train', 'test'])

plt.clf()
x = np.copy(std_data["train"][:,0])
x.sort()
plt.plot(x, my_func(x), marker="o")
plt.savefig("IMG_train_data.pdf")

x = np.copy(std_data["test"][:,0])
x.sort()
plt.plot(x, my_func(x), marker="o")
plt.savefig("IMG_test_data.pdf")

print("Pass here")

