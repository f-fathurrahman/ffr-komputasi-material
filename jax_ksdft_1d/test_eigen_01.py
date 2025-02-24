import jax
import jax.numpy as jnp
from jax.scipy.linalg import qr

jax.config.update("jax_enable_x64", True)

def orthonormalize(x):
    Q, _ = qr(x, mode="economic")
    return Q

# x is assumed to be orthonormalized
def obj_func(x, H):
    Hsub = x.T @ H @ x
    return jnp.trace(Hsub)

# manual gradient
def calc_grad(x, H):
    Hsub = x.T @ H @ x
    g = H @ x - x @ Hsub
    return g


# Optimization step using manual gradient descent
def optimize_eigenvectors(H, x, lr=0.1, steps=500):
    m, n = H.shape
    assert m == n
    
    def loss_fn(x):
        return obj_func(x, H)
    
    loss_grad = jax.grad(loss_fn)
    
    F_old = 0.0
    for iterMin in range(steps):
        g = loss_grad(x) # calc_grad(x, H)
        F = loss_fn(x)
        ΔF = abs(F - F_old)
        print(f"iterMin={iterMin} F={F} ΔF={ΔF}")
        if ΔF < 1e-15:
            print("Converged!")
            break
        #
        x -= lr * g
        x = orthonormalize(x)
        F_old = F
    
    return x



RAND_KEY = jax.random.PRNGKey(1234) # Random seed is explicit in JAX

N = 10
Nstates = 2

H = jax.random.normal(RAND_KEY, shape=(N,N))
H += -10*jnp.eye(N)
H = (H.T + H)/2

# eigenvectors
x = jax.random.normal(RAND_KEY, shape=(N,Nstates))
x = orthonormalize(x)

x_opt = optimize_eigenvectors(H, x)
Hsub = x_opt.T @ H @ x_opt
evals_Hsub = jnp.linalg.eigvalsh(Hsub)
evals_H = jnp.linalg.eigvalsh(H)
print("Some lowest eigenvalues:")
for i in range(Nstates):
    print(f"{i+1}: {evals_Hsub[i]} {evals_H[i]}")

