import numpy as np

import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp


# Define the model
def linear_model(x, w, b):
    return jnp.dot(x, w) + b

def loss(y, labels):
    return jnp.mean((y - labels)**2)

# compute gradients
# w and b are the model parameters
# data is tuple containing input (index-0) and
# targets or labels (index-1)
def loss_wrapper(w, b, data):
    features = data[0]
    labels = data[1]
    y = linear_model(features, w, b)
    return loss(y, labels)

w_true = np.array([0.1])
b_true = 0.3

Ndata = 20
x = np.linspace(-5.0, 5.0, Ndata)
y = np.zeros(Ndata)
for i in range(Ndata):
    y[i] = linear_model(x[i], w_true, b_true)

print("y = ", y)
print("y.shape = ", y.shape)
