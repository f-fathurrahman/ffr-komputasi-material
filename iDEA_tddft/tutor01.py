import numpy as np
import matplotlib.pyplot as plt

import copy
import itertools
import string
import functools
import scipy.sparse as sps
import scipy.sparse.linalg as spsla


import iDEA


#import iDEA.system
#import iDEA.state
#import iDEA.methods.non_interacting

def _estimate_level(s: iDEA.system.System, k: int) -> int:
    r"""
    Estimate the solution to the Schrodinger equation needed to eachive given antisymetric energy state.

    | Args:
    |     s: iDEA.system.System, System object.
    |     k: int, Target energy state.

    | Returns:
    |     level: int, Extimate of level of excitement.
    """
    return (abs(s.up_count - s.down_count) + 1) ** 2 * s.count * (k + 1)


def _permutation_parity(p):
    r"""
    Compute the permulation paritiy of a given permutation.

    | Args:
    |     p: tuple, Permutation.

    | Returns:
    |     parity: float, Permutation parity.
    """
    p = list(p)
    parity = 1
    for i in range(0, len(p) - 1):
        if p[i] != i:
            parity *= -1
            mn = min(range(i, len(p)), key=p.__getitem__)
            p[i], p[mn] = p[mn], p[i]
    return parity



def antisymmetrize(s, spaces, spins, energies):
    r"""
    Antisymmetrize the solution to the Schrodinger equation.

    | Args:
    |     s: iDEA.system.System, System object.
    |     spaces: np.ndarray, Spatial parts of the wavefunction.
    |     spins: np.ndarray, Spin parts of the wavefunction.
    |     energies: np.ndarray, Energies.

    | Returns:
    |     fulls: np.ndarray, Full anantisymmetrized wavefunction.
    |     spaces: np.ndarray, Spatial parts of the wavefunction.
    |     spins: np.ndarray, Spin parts of the wavefunction.
    |     energies: np.ndarray, Energies.

    """
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
    fulls = np.einsum(st, spaces, spins)
    L = list(zip(list(range(0, s.count * 2, 2)), list(range(1, s.count * 2, 2))))
    perms = itertools.permutations(list(range(s.count)))
    fulls_copy = copy.deepcopy(fulls)
    fulls = np.zeros_like(fulls)
    for p in perms:
        indices = list(itertools.chain(*[L[e] for e in p]))
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

    return fulls, spaces, spins, energies




def my_interacting_hamiltonian(s: iDEA.system.System) -> sps.dia_matrix:

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
    terms = generate_terms(sps.kron, h, I, s.count)
    H0 = sps.dia_matrix((s.x.shape[0] ** s.count,) * 2, dtype=float)
    for term in terms:
        H0 += term

    # Add the interaction part of the many-body Hamiltonian
    symbols = string.ascii_lowercase + string.ascii_uppercase
    if s.count > 1:
        indices = ",".join(
            ["".join(c) for c in itertools.combinations(symbols[: s.count], 2)]
        )
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

    return H



def my_interacting_solve(
    s: iDEA.system.System, H: np.ndarray = None, k: int = 0, level=None
) -> iDEA.state.ManyBodyState:
    r"""
    Solves the interacting Schrodinger equation of the given system.

    | Args:
    |     s: iDEA.system.System, System object.
    |     H: np.ndarray, Hamiltonian [If None this will be computed from s]. (default = None)
    |     k: int, Energy state to solve for. (default = 0, the ground-state)
    |     level: int. Max level of excitation to use when solving the Schrodinger equation.

    | Returns:
    |     state: iDEA.state.ManyBodyState, Solved state.
    """
    # Construct the many-body state.
    state = iDEA.state.ManyBodyState()

    # Construct the Hamiltonian.
    if H is None:
        H = my_interacting_hamiltonian(s)

    # Estimate the level of excitation.
    if level is None:
        level = _estimate_level(s, k)

    # Solve the many-body Schrodinger equation.
    print("iDEA.methods.interacting.solve: solving eigenproblem...")
    energies, spaces = spsla.eigsh(H.tocsr(), k=level, which="SA")
    # GPU stuffs are removed

    # Reshape and normalise the solutions.
    spaces = spaces.reshape((s.x.shape[0],) * s.count + (spaces.shape[-1],))
    for j in range(spaces.shape[-1]):
        spaces[..., j] = spaces[..., j] / np.sqrt(
            np.sum(abs(spaces[..., j]) ** 2) * s.dx**s.count
        )

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
    fulls, spaces, spins, energies = antisymmetrize(s, spaces, spins, energies)

    # Populate the state.
    state.space = spaces[..., k]
    state.spin = spins[..., k]
    state.full = fulls[..., k]
    state.energy = energies[k]

    return state



# A predefined system
#atom = iDEA.system.systems.atom

__x2 = np.linspace(-20, 20, 300)
atom = iDEA.system.System(
    __x2, -2.0 / (abs(__x2) + 1.0), iDEA.interactions.softened_interaction(__x2), "ud"
)


ground_state = my_interacting_solve(atom, k=0)