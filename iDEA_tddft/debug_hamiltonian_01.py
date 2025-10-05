import functools
import string
import itertools

import numpy as np
import scipy.sparse as sps

import iDEA

import matplotlib.pyplot as plt

import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")

import matplotlib
matplotlib.style.use("dark_background")

__x2 = np.linspace(-20, 20, 30)
s = iDEA.system.System(
    __x2, -2.0 / (abs(__x2) + 1.0), iDEA.interactions.softened_interaction(__x2), "ud"
)


# Construct the non-interacting part of the many-body Hamiltonian
h = iDEA.methods.non_interacting.hamiltonian(s)[0]

plt.matshow(h);

# Convert to dia matrix.
# XXX: Why not using CSC format?

h = sps.dia_matrix(h)
I = sps.identity(s.x.shape[0], format="dia")

partial_operators = lambda A, B, k, n: (
    A if i + k == n - 1 else B for i in range(n)
)
fold_partial_operators = lambda f, po: functools.reduce(
    lambda acc, val: f(val, acc, format="dia"), po
)
generate_terms = lambda f, A, B, n: (
    fold_partial_operators(f, partial_operators(A, B, k, n)) for k in range(n)
)

# size of H0 is (Npoints^Nelectrons, Npoints^Nelectrons)

H0 = sps.dia_matrix((s.x.shape[0] ** s.count,) * 2, dtype=float)
for term in terms:
    H0 += term

# Add the interaction part of the many-body Hamiltonian
symbols = string.ascii_lowercase + string.ascii_uppercase
if s.count > 1:
    indices = ",".join(
        ["".join(c) for c in itertools.combinations(symbols[: s.count], 2)]
    )
    print("indices = ", indices)
    U = np.log(
        np.einsum(
            indices + "->" + symbols[: s.count],
            *(np.exp(s.v_int),) * int(s.count * (s.count - 1) / 2)
        )
    )
    U = sps.diags(U.reshape((H0.shape[0])), format="dia")
else:
    U = 0.0

for c in itertools.combinations(symbols[: 4], 4):
    print("c = ", c)

# Construct the total many-body Hamiltonian
H = H0 + U

print("type of H = ", type(H))
