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

loss_grad = jax.grad(loss_wrapper, (0, 1))


w_true = np.array([0.1])
b_true = 0.3

Ndata = 20
x = np.linspace(-5.0, 5.0, Ndata)
# Number of feature is 1, we add newaxis to x
x = x[:,np.newaxis]
y = linear_model(x, w_true, b_true)

loss_progress = []
η = 1e-1
data = (x, y)
# Guess value of w and b
w = np.array([1.0])
b = 0.0
for i in range(2000):
    grad = loss_grad(w, b, data)
    w -= η * grad[0]
    b -= η * grad[1]
    L = loss_wrapper(w, b, data)
    print("iteration = {:d} loss = {:f}".format(i, L))
    loss_progress.append(L)
    if L <= 1e-10:
        print("Converged")
        break

print("loss_progress[-1] = ", loss_progress[-1])
print("Final w = ", w)
print("Final b = ", b)

