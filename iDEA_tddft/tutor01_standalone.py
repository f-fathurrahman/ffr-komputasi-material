import numpy as np
import matplotlib.pyplot as plt

import warnings
import copy
import itertools
import string
import functools
import scipy.sparse as sps
import scipy.sparse.linalg as spsla

from abc import ABC as Interface

class State(Interface):
    """Interface class representing a static state."""


class ArrayPlaceholder:
    r"""Array Placeholder."""


class ManyBodyState(State):
    """State of interacting particles."""

    def __init__(
        self, space: np.ndarray = None, spin: np.ndarray = None, full=None, energy=None
    ):
        r"""
        State of particles in a many-body state.

        This is described by a spatial part
        .. math:: \psi(x_1,x_2,\dots,x_N)
        on the spatial grid, and a spin
        part on the spin grid
        .. math:: \chi(\sigma_1,\sigma_2,\dots,\sigma_N).
        These are NOT necessarily antisymmetric states,
        they can be combined using the antisymmetrisation operaration to produce the full
        wavefunction
        .. math:: \Psi(x_1,\sigma_1,x_2,\sigma_2,\dots,x_N,\sigma_N).

        | Args:
        |     space: np.ndarray, Spatial part of the wavefunction on the spatial grid \psi(x_1,x_2,\dots,x_N). (default = None)
        |     spin: np.ndarray, Spin part of the wavefunction on the spin grid \chi(\sigma_1,\sigma_2,\dots,\sigma_N). (default = None)
        |     full: np.ndarray, Total antisymmetrised wavefunction \Psi(x_1,\sigma_1,x_2,\sigma_2,\dots,x_N,\sigma_N). (default = None)
        |     energy: float, Total energy of the state.
        """
        if space is None:
            self.space = ArrayPlaceholder()
        else:
            self.space = space
        if spin is None:
            self.spin = ArrayPlaceholder()
        else:
            self.spin = spin
        if full is None:
            self.full = ArrayPlaceholder()
        else:
            self.full = full
        if energy is None:
            self.energy = float()
        else:
            self.energy = energy


class System:
    r"""Model system, containing all defining properties."""

    def __init__(
        self,
        x: np.ndarray,
        v_ext: np.ndarray,
        v_int: np.ndarray,
        electrons: str,
        stencil: int = 13,
    ):
        r"""
        Model system, containing all defining properties.

        | Args:
        |     x: np.ndarray, Grid of x values in 1D space.
        |     v_ext: np.ndarray, External potential on the grid of x values.
        |     v_int: np.ndarray, Interaction potential on the grid of x values.
        |     electrons: string, Electrons contained in the system.
        |     stencil: int, Stencil to use for derivatives on the grid of x values. (default = 13)

        | Raises:
        |     AssertionError.
        """
        self.__x = x
        self.__dx = self.x[1] - self.x[0]
        self.v_ext = v_ext
        self.v_int = v_int
        self.__electrons = electrons
        self.count = len(electrons)
        self.up_count = electrons.count("u")
        self.down_count = electrons.count("d")
        self.stencil = stencil
        self.check()

    def check(self):
        r"""Performs checks on system properties. Raises AssertionError if any check fails."""
        assert (
            type(self.x) == np.ndarray
        ), f"x grid is not of type np.ndarray, got {type(self.x)} instead."
        assert (
            type(self.v_ext) == np.ndarray
        ), f"v_ext is not of type np.ndarray, got {type(self.v_ext)} instead."
        assert (
            type(self.v_int) == np.ndarray
        ), f"v_int is not of type np.ndarray, got {type(self.v_int)} instead."
        assert (
            type(self.count) == int
        ), f"count is not of type int, got {type(self.NE)} instead."
        assert (
            len(self.x.shape) == 1
        ), f"x grid is not a 1D array, got {len(self.x.shape)}D array instead."
        assert (
            len(self.v_ext.shape) == 1
        ), f"v_ext is not a 1D array, got {len(self.v_ext.shape)}D array instead."
        assert (
            len(self.v_int.shape) == 2
        ), f"v_int is not a 2D array, got {len(self.v_int.shape)}D array instead."
        assert (
            self.x.shape == self.v_ext.shape
        ), f"x grid and v_ext arrays are not the same shape, got x.shape = {self.x.shape} and v_ext.shape = {self.v_ext.shape} instead."
        assert (
            self.x.shape[0] == self.v_int.shape[0]
            and self.x.shape[0] == self.v_int.shape[1]
        ), "v_int is not of the correct shape, got shape {self.v_int.shape} instead."
        assert self.count >= 0, f"count is not positive."
        assert set(self.electrons).issubset(
            set(["u", "d"])
        ), f"Electrons must have only up or down spin, e.g 'uudd'. Got {self.electrons} instead"
        assert (
            self.count == self.up_count + self.down_count
        ), f"Electrons must obay up_count + down_count = count."
        assert self.stencil in [
            3,
            5,
            7,
            9,
            11,
            13,
        ], f"stencil must be one of [3,5,7,9,11,13], got {self.stencil} instead."

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, value):
        self.__x = value
        self.__dx = self.__x[1] - self.__x[0]
        warnings.warn(
            "x grid has been changed: dx has been recomputed, please update v_ext and v_int on this grid."
        )

    @x.deleter
    def x(self):
        del self.__x

    @property
    def dx(self):
        return self.__dx

    @dx.setter
    def dx(self, value):
        raise AttributeError(
            "cannot set dx directly: set the x grid and dx will be updated automatically."
        )

    @dx.deleter
    def dx(self):
        del self.__dx

    @property
    def electrons(self):
        return self.__electrons

    @electrons.setter
    def electrons(self, value):
        self.__electrons = value
        self.count = len(value)
        self.up_count = value.count("u")
        self.down_count = value.count("d")

    @electrons.deleter
    def electrons(self):
        del self.__electrons

    def __str__(self):
        return f"System: x = np.array([{self.x[0]:.3f},...,{self.x[-1]:.3f}]), dx = {self.dx:.4f}..., v_ext = np.array([{self.v_ext[0]:.3f},...,{self.v_ext[-1]:.3f}]), electrons = {self.electrons}"


# XXX: rename to one-electron Hamiltonian?
def non_interacting_hamiltonian(
    s: System,
    up_n: np.ndarray = None,
    down_n: np.ndarray = None,
    up_p: np.ndarray = None,
    down_p: np.ndarray = None,
    K: np.ndarray = None,
    Vext: np.ndarray = None,
) -> np.ndarray:
    r"""
    Compute the Hamiltonian from the kinetic and potential terms.

    | Args:
    |     s: System, System object.
    |     up_n: np.ndarray, Charge density of up electrons.
    |     down_n: np.ndarray, Charge density of down electrons.
    |     up_p: np.ndarray, Charge density matrix of up electrons.
    |     down_p: np.ndarray, Charge density matrix of down electrons.
    |     K: np.ndarray, Single-particle kinetic energy operator [If None this will be computed from s]. (default = None)
    |     Vext: np.ndarray, Potential energy operator [If None this will be computed from s]. (default = None)

    | Returns:
    |     H: np.ndarray, Hamiltonian, up Hamiltonian, down Hamiltonian.
    """
    if K is None:
        K = kinetic_energy_operator(s)
    if Vext is None:
        Vext = external_potential_operator(s)
    H = K + Vext
    return H, H, H



def kinetic_energy_operator(s: System) -> np.ndarray:
    r"""
    Compute single-particle kinetic energy operator as a matrix.

    This is built using a given number of finite differences to represent the second derivative.
    The number of differences taken is defined in s.stencil.

    | Args:
    |     s: System, System object.

    | Returns:
    |     K: np.ndarray, Kintetic energy operator.
    """
    if s.stencil == 3:
        sd = 1.0 * np.array([1, -2, 1], dtype=float) / s.dx**2
        sdi = (-1, 0, 1)
    elif s.stencil == 5:
        sd = 1.0 / 12.0 * np.array([-1, 16, -30, 16, -1], dtype=float) / s.dx**2
        sdi = (-2, -1, 0, 1, 2)
    elif s.stencil == 7:
        sd = (
            1.0
            / 180.0
            * np.array([2, -27, 270, -490, 270, -27, 2], dtype=float)
            / s.dx**2
        )
        sdi = (-3, -2, -1, 0, 1, 2, 3)
    elif s.stencil == 9:
        sd = (
            1.0
            / 5040.0
            * np.array(
                [-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9], dtype=float
            )
            / s.dx**2
        )
        sdi = (-4, -3, -2, -1, 0, 1, 2, 3, 4)
    elif s.stencil == 11:
        sd = (
            1.0
            / 25200.0
            * np.array(
                [8, -125, 1000, -6000, 42000, -73766, 42000, -6000, 1000, -125, 8],
                dtype=float,
            )
            / s.dx**2
        )
        sdi = (-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5)
    elif s.stencil == 13:
        sd = (
            1.0
            / 831600.0
            * np.array(
                [
                    -50,
                    864,
                    -7425,
                    44000,
                    -222750,
                    1425600,
                    -2480478,
                    1425600,
                    -222750,
                    44000,
                    -7425,
                    864,
                    -50,
                ],
                dtype=float,
            )
            / s.dx**2
        )
        sdi = (-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6)
    second_derivative = np.zeros((s.x.shape[0], s.x.shape[0]))
    for i in range(len(sdi)):
        second_derivative += np.diag(
            np.full(
                np.diag(np.zeros((s.x.shape[0], s.x.shape[0])), k=sdi[i]).shape[0],
                sd[i],
            ),
            k=sdi[i],
        )
    K = -0.5 * second_derivative
    return K


def external_potential_operator(s: System) -> np.ndarray:
    Vext = np.diag(s.v_ext)
    return Vext




def _estimate_level(s: System, k: int) -> int:
    r"""
    Estimate the solution to the Schrodinger equation needed to eachive given antisymetric energy state.

    | Args:
    |     s: System, System object.
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
    |     s: System, System object.
    |     spaces: np.ndarray, Spatial parts of the wavefunction.
    |     spins: np.ndarray, Spin parts of the wavefunction.
    |     energies: np.ndarray, Energies.

    | Returns:
    |     fulls: np.ndarray, Full antisymmetrized wavefunction.
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




def my_interacting_hamiltonian(s: System) -> sps.dia_matrix:

    # Construct the non-interacting part of the many-body Hamiltonian
    h = non_interacting_hamiltonian(s)[0] # only use the first returned result, total Hamiltonian
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
    s: System, H: np.ndarray = None, k: int = 0, level=None
) -> ManyBodyState:
    r"""
    Solves the interacting Schrodinger equation of the given system.

    | Args:
    |     s: System, System object.
    |     H: np.ndarray, Hamiltonian [If None this will be computed from s]. (default = None)
    |     k: int, Energy state to solve for. (default = 0, the ground-state)
    |     level: int. Max level of excitation to use when solving the Schrodinger equation.

    | Returns:
    |     state: ManyBodyState, Solved state.
    """
    # Construct the many-body state.
    state = ManyBodyState()

    # Construct the Hamiltonian.
    if H is None:
        H = my_interacting_hamiltonian(s)
    print("Hamiltonian size = ", H.shape)

    # Estimate the level of excitation.
    if level is None:
        level = _estimate_level(s, k)
    print("level = ", level)

    # Solve the many-body Schrodinger equation.
    print("my_interacting_solve: solving eigenproblem...")
    energies, spaces = spsla.eigsh(H.tocsr(), k=level, which="SA")
    print("energies raw = ", energies)

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




def softened_interaction(
    x: np.ndarray, strength: float = 1.0, softening: float = 1.0
) -> np.ndarray:
    v_int = np.zeros((x.shape[0], x.shape[0]), dtype="float")
    for i in range(x.shape[0]):
        for j in range(x.shape[0]):
            v_int[i, j] = strength / (abs(x[i] - x[j]) + softening)
    return v_int


def softened_interaction_alternative(
    x: np.ndarray, strength: float = 1.0, softening: float = 1.0
) -> np.ndarray:
    v_int = np.zeros((x.shape[0], x.shape[0]), dtype="float")
    for i in range(x.shape[0]):
        for j in range(x.shape[0]):
            v_int[i, j] = strength / np.sqrt(((x[i] - x[j]) ** 2 + softening))
    return v_int


#__x2 = np.linspace(-20, 20, 50)
#atom = System(
#    __x2, -2.0 / (abs(__x2) + 1.0), softened_interaction(__x2), "ud"
#)

#ω = 1.0
#x = np.linspace(-10, 10, 150)
#v_ext = 0.5 * ω**2 * x**2
#v_int = softened_interaction_alternative(x)
#atom = System(x, v_ext, v_int, electrons='ud') # also test ud
#
#ground_state = my_interacting_solve(atom, k=0)
#print("Total energy = ", ground_state.energy)



ω = 1.0
x = np.linspace(-8, 8, 20)
v_ext = 0.5 * ω**2 * x**2
v_int = softened_interaction_alternative(x)
atom = System(x, v_ext, v_int, electrons='udu', stencil=3)

ground_state = my_interacting_solve(atom, k=0)
print("Total energy = ", ground_state.energy)


