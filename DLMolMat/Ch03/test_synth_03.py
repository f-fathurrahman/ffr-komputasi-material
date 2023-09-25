import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("seaborn-v0_8-darkgrid")
plt.rcParams["savefig.dpi"] = 150

import jax
import jax.numpy as jnp

# define our loss function
@jax.jit
def loss(w, x, y):
    return jnp.mean((y - jnp.dot(x, w)) ** 2)

np.random.seed(1234)

# generate data from polynomial
N = 20
syn_x = np.linspace(-3, 3, N)
# create feature matrix or design matrix
syn_features = np.vstack([syn_x**3, syn_x**2, syn_x, np.ones_like(syn_x)]).T
# sequence: x^3, x^2, x^1, and x^0

# The labels (function values)
syn_labels = 2*syn_x**3 - syn_x**2 + syn_x - 1
# coeffs w

#
# split data into train/test
#
# Train data: 1/4 first and 1/4 last
indices = list(range(0, N // 4)) + list(range(3 * N // 4, N)) 
# test data (between 1/4 and 3/4)
test_indices = list(range(N // 4, 3 * N // 4))


# Let's see an example where the feature number is the same but they aren't
# perfectly correlated with labels, meaning we cannot match the labels even
# if there was no noise.
syn_features = np.vstack(
    [syn_x**2, syn_x, np.exp(-(syn_x**2)), np.cos(syn_x), np.ones_like(syn_x)]
).T

train_x = syn_features[indices]
train_y = syn_labels[indices]

test_x = syn_features[test_indices]
test_y = syn_labels[test_indices]

# Add noise
train_y = train_y + np.random.normal(scale=5, size=train_y.shape)

w, *_ = np.linalg.lstsq(train_x, train_y, rcond=-1)
print("w = ", w)

plt.clf()
plt.plot(syn_x[indices], train_y, "o", label="Training labels")
plt.plot(syn_x[test_indices], test_y, "o", label="Testing labels")
plt.ylim(-80, 80)
plt.plot(syn_x, jnp.dot(syn_features, w), label="Fit Model")
plt.plot(syn_x, syn_labels, "--", label="Ground Truth")
plt.text(0, -20, f"Training Loss {loss(w, train_x, train_y):.2f}")
plt.text(0, -30, f"Testing Loss {loss(w, test_x, test_y):.2f}")
plt.legend(loc="upper left")
plt.title("Noise, Imperfectly Correlated Features")
plt.savefig("IMG_synth_noise_04.png")
plt.show()