import numpy as np
import jax
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

# loss grad is a new function
loss_grad = jax.grad(loss_wrapper, (0, 1))

w_true = np.array([0.1, 0.3, 0.4])
b_true = 0.1

# Some data
x = np.array([1.0, 0.1, 2.5]) # we only have one instance

# The target, generated using true model parameter
y_true = w_true[0]*x[0] + w_true[1]*x[1] + w_true[2]*x[2] + b_true
print("y_true = ", y_true)

# Guess the value of parameters
w = np.array([0.1, 0.3, 0.4])
b = 0.1

y_pred = linear_model(x, w, b)
print("y_pred = ", y_pred)

print("loss = ", loss_wrapper(w, b, (x, y_true)))

dLdx, dLdb = loss_grad(w, b, (x, y_true))
print("loss grad dLdx = ", dLdx)
print("loss grad dLdb = ", dLdb)
