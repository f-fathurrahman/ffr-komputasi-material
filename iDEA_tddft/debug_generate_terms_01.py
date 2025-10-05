import functools
import string
import itertools

import numpy as np
import scipy.sparse as sps

import iDEA

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("dark_background")

Npoints = 5
__x2 = np.linspace(-20, 20, Npoints)
s = iDEA.system.System(
    __x2, -2.0 / (abs(__x2) + 1.0), iDEA.interactions.softened_interaction(__x2), "udu", stencil=3
)

print("Number of electrons = ", s.count)
print("electrons config = ", s.electrons)
print("No. of spin up = ", s.up_count)
print("No. of spin down = ", s.down_count)


# Construct the non-interacting part of the many-body Hamiltonian
h = iDEA.methods.non_interacting.hamiltonian(s)[0]


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

def func_partial_operators(A, B, k, n):
    ret_vals = []
    for i in range(n):
        if i + k == n - 1:
            ret_vals.append("h")
        else:
            ret_vals.append("I")
    return ret_vals

Ntry = s.count
terms_str = [ func_partial_operators(h, I, k, Ntry) for k in range(Ntry) ]
print("terms_str = ", terms_str)

terms = generate_terms(sps.kron, h, I, s.count)

#TERMS = []
#for term in terms:
#    TERMS.append(term)

H0 = sps.dia_matrix((s.x.shape[0] ** s.count,) * 2, dtype=float)
for term in terms:
    H0 += term

# size of H0 is (Npoints^Nelectrons, Npoints^Nelectrons)

# Manual construction
#HH = sps.kron(I, h) + sps.kron(h, I) # two electrons
HH = sps.kron( sps.kron(I, I), h ) + \
     sps.kron( sps.kron(I, h), I ) + \
     sps.kron( sps.kron(h, I), I )
