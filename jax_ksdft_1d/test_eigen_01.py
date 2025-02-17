import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

RAND_KEY = jax.random.PRNGKey(1234) # Random seed is explicit in JAX

N = 5
Nstates = 2

H = jax.random.normal(RAND_KEY, shape=(N,N))
H += -10*jnp.eye(N)

# eigenvectors
x = jax.random.normal(RAND_KEY, shape=(N,Nstates))

def orthonormalize(x):
    Q, _ = qr(x)
    return Q

def obj_func(x, H):
    Hsub = x.T @ H @ x
    return jnp.trace(Hsub)
