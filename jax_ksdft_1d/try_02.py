import jax
import jax.numpy as jnp
from jax.scipy.linalg import qr
from jax import grad, jit

jax.config.update("jax_enable_x64", True)

# Function to compute the Rayleigh-Ritz quotient for k eigenvectors
def rayleigh_quotient(V, A):
    return jnp.trace(V.T @ A @ V)  # Sum of k Rayleigh quotients

# Ensure V is orthonormal using QR decomposition
def project_onto_stiefel(V):
    Q, _ = qr(V, mode="economic")
    return Q

# Optimization step using manual gradient descent
def optimize_eigenvectors(A, k, lr=0.1, steps=500):
    m, n = A.shape
    key = jax.random.PRNGKey(42)
    V = jax.random.normal(key, (n, k))
    V = project_onto_stiefel(V)
    
    @jit
    def loss_fn(V):
        return rayleigh_quotient(V, A)
    
    loss_grad = jit(grad(loss_fn))
    
    for _ in range(steps):
        grads = loss_grad(V)
        #
        norm_grad = 0.0
        for ist in range(n):
            norm_grad += jnp.dot(grads[:,ist], grads[:,ist])
        print("norm grad = ", norm_grad)
        #
        V -= lr * grads  # Gradient descent step
        V = project_onto_stiefel(V)  # Re-orthonormalize
    
    return V

# Example usage
# Example usage
N = 5
Nstates = 2

RAND_KEY = jax.random.PRNGKey(1234) # Random seed is explicit in JAX

H = jax.random.normal(RAND_KEY, shape=(N,N))
H += -10*jnp.eye(N) # make it diagonally dominant
H = (H + H.T)/2 # make it symmetric

V_opt = optimize_eigenvectors(H, Nstates)
Hsub = V_opt.T @ H @ V_opt

evals_Hsub = jnp.linalg.eigvalsh(Hsub)
evals_H = jnp.linalg.eigvalsh(H)
print("Some lowest eigenvalues:")
for i in range(Nstates):
    print(f"{i+1}: {evals_Hsub[i]} {evals_H[i]}")

