# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %% [markdown]
# # Tensor Manipulation: Psi4 and NumPy manipulation routines
# Contracting tensors together forms the core of the Psi4NumPy project. First let us consider the popluar [Einstein Summation Notation](https://en.wikipedia.org/wiki/Einstein_notation) which allows for very succinct descriptions of a given tensor contraction.
#
# For example, let us consider a [inner (dot) product](https://en.wikipedia.org/wiki/Dot_product):
# $$c = \sum_{ij} A_{ij} * B_{ij}$$
#
# With the Einstein convention, all indices that are repeated are considered summed over, and the explicit summation symbol is dropped:
# $$c = A_{ij} * B_{ij}$$
#
# This can be extended to [matrix multiplication](https://en.wikipedia.org/wiki/Matrix_multiplication):
# \begin{align}
# \rm{Conventional}\;\;\;  C_{ik} &= \sum_{j} A_{ij} * B_{jk} \\
# \rm{Einstein}\;\;\;  C &= A_{ij} * B_{jk} \\
# \end{align}
#
# Where the $C$ matrix has *implied* indices of $C_{ik}$ as the only repeated index is $j$.
#
# However, there are many cases where this notation fails. Thus we often use the generalized Einstein convention. To demonstrate let us examine a [Hadamard product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)):
# $$C_{ij} = \sum_{ij} A_{ij} * B_{ij}$$
#
#
# This operation is nearly identical to the dot product above, and is not able to be written in pure Einstein convention. The generalized convention allows for the use of indices on the left hand side of the equation:
# $$C_{ij} = A_{ij} * B_{ij}$$
#
# Usually it should be apparent within the context the exact meaning of a given expression.
#
# Finally we also make use of Matrix notation:
# \begin{align}
# {\rm Matrix}\;\;\;  \bf{D} &= \bf{A B C} \\
# {\rm Einstein}\;\;\;  D_{il} &= A_{ij} B_{jk} C_{kl}
# \end{align}
#
# Note that this notation is signified by the use of bold characters to denote matrices and consecutive matrices next to each other imply a chain of matrix multiplications! 

# %% [markdown]
# ## Einsum
#
# To perform most operations we turn to [NumPy's einsum function](https://docs.scipy.org/doc/numpy/reference/generated/numpy.einsum.html) which allows the Einsten convention as an input. In addition to being much easier to read, manipulate, and change, it is also much more efficient that a pure Python implementation.
#
# To begin let us consider the construction of the following tensor (which you may recognize):
# $$G_{pq} = 2.0 * I_{pqrs} D_{rs} - 1.0 * I_{prqs} D_{rs}$$ 
#
# First let us import our normal suite of modules:

# %%
import numpy as np
import psi4
import time

# %% [markdown]
# We can then use conventional Python loops and einsum to perform the same task. Keep size relatively small as these 4-index tensors grow very quickly in size.

# %%
size = 20

if size > 30:
    raise Exception("Size must be smaller than 30.")
D = np.random.rand(size, size)
I = np.random.rand(size, size, size, size)

# Build the fock matrix using loops, while keeping track of time
tstart_loop = time.time()
Gloop = np.zeros((size, size))
for p in range(size):
    for q in range(size):
        for r in range(size):
            for s in range(size):
                Gloop[p, q] += 2 * I[p, q, r, s] * D[r, s]
                Gloop[p, q] -=     I[p, r, q, s] * D[r, s]

g_loop_time = time.time() - tstart_loop

# Build the fock matrix using einsum, while keeping track of time
tstart_einsum = time.time()
J = np.einsum('pqrs,rs', I, D, optimize=True)
K = np.einsum('prqs,rs', I, D, optimize=True)
G = 2 * J - K

einsum_time = time.time() - tstart_einsum

# Make sure the correct answer is obtained
print('The loop and einsum fock builds match:    %s\n' % np.allclose(G, Gloop))
# Print out relative times for explicit loop vs einsum Fock builds
print('Time for loop G build:   %14.4f seconds' % g_loop_time)
print('Time for einsum G build: %14.4f seconds' % einsum_time)
print('G builds with einsum are {:3.4f} times faster than Python loops!'.format(g_loop_time / einsum_time))

# %% [markdown]
# As you can see, the einsum function is considerably faster than the pure Python loops and, in this author's opinion, much cleaner and easier to use.

# %% [markdown]
# ## Dot
#
# Now let us turn our attention to a more canonical matrix multiplication example such as:
# $$D_{il} = A_{ij} B_{jk} C_{kl}$$
#
# We could perform this operation using einsum; however, matrix multiplication is an extremely common operation in all branches of linear algebra. Thus, these functions have been optimized to be more efficient than the `einsum` function. The matrix product will explicitly compute the following operation:
# $$C_{ij} = A_{ij} * B_{ij}$$
#
# This can be called with [NumPy's dot function](https://docs.scipy.org/doc/numpy/reference/generated/numpy.dot.html#numpy.dot).

# %%
size = 200
A = np.random.rand(size, size)
B = np.random.rand(size, size)
C = np.random.rand(size, size)

# First compute the pair product
tmp_dot = np.dot(A, B)
tmp_einsum = np.einsum('ij,jk->ik', A, B, optimize=True)
print("Pair product allclose: %s" % np.allclose(tmp_dot, tmp_einsum))

# %% [markdown]
# Now that we have proved exactly what the dot product does, let us consider the full chain and do a timing comparison:

# %%
D_dot = np.dot(A, B).dot(C)
D_einsum = np.einsum('ij,jk,kl->il', A, B, C, optimize=True)
print("Chain multiplication allclose: %s" % np.allclose(D_dot, D_einsum))

print("\nnp.dot time:")
# %timeit np.dot(A, B).dot(C)

print("\nnp.einsum time")
# no optimization here for illustrative purposes!
# %timeit np.einsum('ij,jk,kl->il', A, B, C)

# %% [markdown]
# On most machines the `np.dot` times are roughly ~2,000 times faster. The reason is twofold:
#  - The `np.dot` routines typically call [Basic Linear Algebra Subprograms (BLAS)](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms). The BLAS routines are highly optimized and threaded versions of the code.
#  - The `np.einsum` code will not factorize the operation by default; Thus, the overall cost is ${\cal O}(N^4)$ (as there are four indices) rather than the factored $(\bf{A B}) \bf{C}$ which runs ${\cal O}(N^3)$.
#  
# The first issue is difficult to overcome; however, the second issue can be resolved by the following:

# %%
print("np.einsum factorized time:")
# no optimization here for illustrative purposes!
# %timeit np.einsum('ik,kl->il', np.einsum('ij,jk->ik', A, B), C)

# %% [markdown]
# On most machines the factorized `einsum` expression is only ~10 times slower than `np.dot`. While a massive improvement, this is a clear demonstration the BLAS usage is usually recommended. It is a tradeoff between speed and readability. The Psi4NumPy project tends to lean toward `einsum` usage except in case where the benefit is too large to pass up.
#
# Starting in NumPy 1.12, the [einsum function](https://docs.scipy.org/doc/numpy/reference/generated/numpy.einsum.html) has a `optimize` flag which will automatically factorize the einsum code for you using a greedy algorithm, leading to considerable speedups at almost no cost:

# %%
print("\nnp.einsum optimized time")
# %timeit np.einsum('ij,jk,kl->il', A, B, C, optimize=True)

# %% [markdown]
# In this example, using `optimize=True` for automatic factorization is only 25% slower than `np.dot`. Furthermore, it is ~5 times faster than factorizing the expression by hand, which represents a very good trade-off between speed and readability. When unsure, `optimize=True` is strongly recommended.

# %% [markdown]
# ## Complex tensor manipulations
# Let us consider a popular index transformation example:
# $$M_{pqrs} = C_{pi} C_{qj} I_{ijkl} C_{rk} C_{sl}$$
#
# Here, a naive `einsum` call would scale like $\mathcal{O}(N^8)$ which translates to an extremely costly computation for all but the smallest $N$.

# %%
# Grab orbitals
size = 15
if size > 15:
    raise Exception("Size must be smaller than 15.")
    
C = np.random.rand(size, size)
I = np.random.rand(size, size, size, size)

# Numpy einsum N^8 transformation.
print("\nStarting Numpy's N^8 transformation...")
n8_tstart = time.time()
# no optimization here for illustrative purposes!
MO_n8 = np.einsum('pI,qJ,pqrs,rK,sL->IJKL', C, C, I, C, C)
n8_time = time.time() - n8_tstart
print("...transformation complete in %.3f seconds." % (n8_time))

# Numpy einsum N^5 transformation.
print("\n\nStarting Numpy's N^5 transformation with einsum...")
n5_tstart = time.time()
# no optimization here for illustrative purposes!
MO_n5 = np.einsum('pA,pqrs->Aqrs', C, I)
MO_n5 = np.einsum('qB,Aqrs->ABrs', C, MO_n5)
MO_n5 = np.einsum('rC,ABrs->ABCs', C, MO_n5)
MO_n5 = np.einsum('sD,ABCs->ABCD', C, MO_n5)
n5_time = time.time() - n5_tstart
print("...transformation complete in %.3f seconds." % n5_time)
print("\nN^5 %4.2f faster than N^8 algorithm!" % (n8_time / n5_time))
print("Allclose: %s" % np.allclose(MO_n8, MO_n5))

# Numpy einsum optimized transformation.
print("\nNow Numpy's optimized transformation...")
n8_tstart = time.time()
MO_n8 = np.einsum('pI,qJ,pqrs,rK,sL->IJKL', C, C, I, C, C, optimize=True)
n8_time_opt = time.time() - n8_tstart
print("...optimized transformation complete in %.3f seconds." % (n8_time_opt))

# Numpy GEMM N^5 transformation.
# Try to figure this one out!
print("\n\nStarting Numpy's N^5 transformation with dot...")
dgemm_tstart = time.time()
MO = np.dot(C.T, I.reshape(size, -1))
MO = np.dot(MO.reshape(-1, size), C)
MO = MO.reshape(size, size, size, size).transpose(1, 0, 3, 2)

MO = np.dot(C.T, MO.reshape(size, -1))
MO = np.dot(MO.reshape(-1, size), C)
MO = MO.reshape(size, size, size, size).transpose(1, 0, 3, 2)
dgemm_time = time.time() - dgemm_tstart
print("...transformation complete in %.3f seconds." % dgemm_time)
print("\nAllclose: %s" % np.allclose(MO_n8, MO))
print("N^5 %4.2f faster than N^8 algorithm!" % (n8_time / dgemm_time))
