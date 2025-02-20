import jax
import jax.numpy as jnp
from jax.scipy.linalg import qr

jax.config.update("jax_enable_x64", True)

N = 5
Nstates = 2

RAND_KEY = jax.random.PRNGKey(1234)
x = jax.random.normal(RAND_KEY, shape=(N,Nstates))

print("Initial: ")
print(x.T @ x)

x, _ = qr(x, mode="economic")
print("After: ")
print(x.T @ x)
