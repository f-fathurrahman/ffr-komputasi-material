import autograd.numpy as anp
import pymanopt

anp.random.seed(42)

dim = 10
manifold = pymanopt.manifolds.Sphere(dim)
matrix = anp.random.normal(size=(dim,dim))
matrix = 0.5*(matrix + matrix.T) # symmetrize

@pymanopt.function.autograd(manifold)
def cost(point):
    # Need negative sign because the algorithm is minimization
    return -point @ matrix @ point

problem = pymanopt.Problem(manifold, cost)
#optimizer = pymanopt.optimizers.SteepestDescent()
optimizer = pymanopt.optimizers.ConjugateGradient()
result = optimizer.run(problem)

eigenvalues, eigenvectors = anp.linalg.eig(matrix)
dominant_eigenvector = eigenvectors[:, eigenvalues.argmax()]

# Evalute the actual largest eigenvalue.
# This is the same as cost function without negative sign.
X = result.point
λ = X.T @ (matrix @ X)
# Actually the same as X @ matrix @ X

