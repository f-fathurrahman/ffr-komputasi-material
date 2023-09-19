import numpy as np
import jax
import jax.numpy as jnp

# Define the model
def linear_model(x, w, b):
    return jnp.dot(x, w) + b

# Some test
x = np.array([1, 0, 2.5])
w = np.array([0.2, -0.5, 0.4])
b = 4.3

print(linear_model(x, w, b))
