import functools
import string
import itertools
import copy

import numpy as np
import scipy.sparse as sps

import iDEA

import matplotlib.pyplot as plt

import matplotlib_inline
matplotlib_inline.backend_inline.set_matplotlib_formats("svg")

import matplotlib
matplotlib.style.use("dark_background")

def _permutation_parity(p):
    p = list(p)
    parity = 1
    for i in range(0, len(p) - 1):
        if p[i] != i:
            parity *= -1
            mn = min(range(i, len(p)), key=p.__getitem__)
            p[i], p[mn] = p[mn], p[i]
    return parity


Npoints = 5
__x2 = np.linspace(-20, 20, Npoints)
s = iDEA.system.System(
    __x2, -2.0 / (abs(__x2) + 1.0), iDEA.interactions.softened_interaction(__x2), "udu", stencil=3
)

# Construct the non-interacting part of the many-body Hamiltonian
h = iDEA.methods.non_interacting.hamiltonian(s)[0]


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
terms = generate_terms(sps.kron, h, I, s.count)

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

# Construct the total many-body Hamiltonian
H = H0 + U

print("type of H = ", type(H))

H0dense = H0.toarray()
Udense = U.toarray()

k = 0 # ground state?
level = (abs(s.up_count - s.down_count) + 1) ** 2 * s.count * (k + 1)
print("level = ", level)


import scipy.sparse.linalg as spsla
energies, spaces = spsla.eigsh(H.tocsr(), k=level, which="SA")

print("energies = ", energies)


# Reshape and normalise the solutions.
print("BEFORE: Shape of spaces (spatial wavefunction) = ", spaces.shape)
spaces = spaces.reshape((s.x.shape[0],) * s.count + (spaces.shape[-1],))
print("AFTER: Shape of spaces (spatial wavefunction) = ", spaces.shape)

for j in range(spaces.shape[-1]):
    cc = np.sum(abs(spaces[..., j])**2) * s.dx**s.count
    spaces[..., j] = spaces[..., j] / np.sqrt(cc)


# Construct the spin part.
symbols = string.ascii_lowercase + string.ascii_uppercase
u = np.array([1, 0])
d = np.array([0, 1])
spin_state = tuple([u if spin == "u" else d for spin in s.electrons])
spin = np.einsum(
    ",".join(symbols[: s.count]) + "->" + "".join(symbols[: s.count]), *spin_state
)
spins = np.zeros(shape=((2,) * s.count + (spaces.shape[-1],)))
for i in range(spaces.shape[-1]):
    spins[..., i] = spin

# Antisymmetrize.
#fulls, spaces, spins, energies = antisymmetrize(s, spaces, spins, energies)


# Perform antisymmetrization.
l = string.ascii_lowercase[: s.count]
L = string.ascii_uppercase[: s.count]
st = (
    l
    + "Y,"
    + L
    + "Y->"
    + "".join([i for sub in list(zip(l, L)) for i in sub])
    + "Y"
)
print("st = ", st)
fulls = np.einsum(st, spaces, spins)
L = list(zip(list(range(0, s.count * 2, 2)), list(range(1, s.count * 2, 2))))
perms = itertools.permutations(list(range(s.count)))
fulls_copy = copy.deepcopy(fulls)
fulls = np.zeros_like(fulls)
for p in perms:
    indices = list(itertools.chain(*[L[e] for e in p]))
    print("permutation = ", p, " parity = ", _permutation_parity(p))
    fulls += _permutation_parity(p) * np.moveaxis(
        fulls_copy, list(range(s.count * 2)), indices
    )

# Filter out zeros.
allowed_fulls = []
allowed_energies = []
allowed_spaces = []
allowed_spins = []
for n in range(fulls.shape[-1]):
    if np.allclose(fulls[..., n], np.zeros(fulls.shape[:-1])):
        pass
    else:
        allowed_fulls.append(fulls[..., n])
        allowed_energies.append(energies[n])
        allowed_spaces.append(spaces[..., n])
        allowed_spins.append(spins[..., n])
fulls = np.moveaxis(np.array(allowed_fulls), 0, -1)
spaces = np.moveaxis(np.array(allowed_spaces), 0, -1)
spins = np.moveaxis(np.array(allowed_spins), 0, -1)
energies = np.array(allowed_energies)

# Normalise.
for k in range(fulls.shape[-1]):
    fulls[..., k] = fulls[..., k] / np.sqrt(
        np.sum(abs(fulls[..., k]) ** 2) * s.dx**s.count
    )

# Filter out duplicates.
allowed_fulls = []
allowed_energies = []
for n in range(fulls.shape[-1] - 1):
    if np.allclose(abs(fulls[..., n]), abs(fulls[..., n + 1])):
        pass
    else:
        allowed_fulls.append(fulls[..., n])
        allowed_energies.append(energies[n])
allowed_fulls.append(fulls[..., -1])
allowed_energies.append(energies[-1])
fulls = np.moveaxis(np.array(allowed_fulls), 0, -1)
spaces = spaces[..., : fulls.shape[-1]]
spins = spins[..., : fulls.shape[-1]]
energies = np.array(allowed_energies)

#return fulls, spaces, spins, energies