import ase
import matplotlib.pyplot as plt

###############################################################################
# Simple example using my_boltztrap2 as a library to calculate the Onsager
# coefficients for a simple cubic structure containing a single, parabolic,
# isotropic band.
# Strict parabolicity and isotropy are broken by the cubic lattice, but hold
# for energies close to the botton of the band.
###############################################################################
import numpy as np
import scipy as sp
import scipy.linalg as la
import spglib

import my_boltztrap2.bandlib
import my_boltztrap2.fite
import my_boltztrap2.sphere
from my_boltztrap2.units import *

# Number of k points along each direction in the input grid
NK = 25
# Effective mass of the only band
EFFM = 1.0
# Minimum energy of the band
OFFSET = 0.5
# Amplification factor for the interpolation
FACTOR = 5

###############################################################################
# Step 1: Create the structure and the band
###############################################################################

# Create a hypotetical monoatomic simple cubic structure with a lattice
# parameter of 5 A.
atoms = ase.Atoms("Si", cell=5 * np.eye(3), pbc=True)
lattvec = atoms.get_cell().T * Angstrom
rlattvec = 2.0 * np.pi * la.inv(lattvec).T
# Create a set of irreducible k points based on that structure
fromspglib = spglib.get_ir_reciprocal_mesh(
    [NK, NK, NK],
    (
        atoms.get_cell(),
        atoms.get_scaled_positions(),
        atoms.get_atomic_numbers(),
    ),
)
indices = np.unique(fromspglib[0]).tolist()
kpoints = fromspglib[1].T / float(NK)
kpoints = kpoints[:, indices]
# Compute the band energies in a purely parabolic model
cartesian = rlattvec @ kpoints
k2 = (cartesian**2).sum(axis=0)
eband = OFFSET + k2 / 2.0 / EFFM


# Weight the band by a bump function to make its derivative zero at the
# boundary of the BZ.
def create_bump(a, b):
    """Return a bump function f(x) equal to zero for |x| > b, equal to one for
    |x| < a, and providing a smooth transition in between.
    """
    if a <= 0.0 or b <= 0.0 or a >= b:
        raise ValueError("a and b must be positive numbers, with b > a")

    # Follow the prescription given by Loring W. Tu in
    # "An Introduction to Manifolds", 2nd Edition, Springer
    def f(t):
        """Auxiliary function used as a building block of the bump function."""
        if t <= 0.0:
            return 0.0
        else:
            return np.exp(-1.0 / t)

    f = np.vectorize(f)

    def nruter(x):
        """One-dimensional bump function."""
        arg = (x * x - a * a) / (b * b - a * a)
        return 1.0 - f(arg) / (f(arg) + f(1.0 - arg))

    return nruter


rmax = rlattvec[0, 0] / 2.0
rmin = 0.8 * rmax
bound = eband[k2 <= rmax * rmax].max()
bump = create_bump(rmin, rmax)
eband = np.minimum(eband, bound)
weight = bump(cartesian[0, :]) * bump(cartesian[1, :]) * bump(cartesian[2, :])
eband = eband * weight + bound * (1.0 - weight)


# Create a minimalistic mock DFTData object from these elements
class MockDFTData:
    """Mock DFTData emulation class for the parabolic band example."""

    def __init__(self):
        """Create a MockDFTData object based on global variables."""
        self.kpoints = np.copy(kpoints.T)
        self.ebands = np.copy(eband.reshape((1, eband.size)))
        self.mommat = None

    def get_lattvec(self):
        """Return the matrix of lattice vectors."""
        return lattvec


data = MockDFTData()

if __name__ == "__main__":
    ###############################################################################
    # Step 2: Obtain the interpolation coefficients
    ###############################################################################
    equivalences = my_boltztrap2.sphere.get_equivalences(
        atoms, None, FACTOR * kpoints.shape[1]
    )
    coeffs = my_boltztrap2.fite.fitde3D(data, equivalences)

    ###############################################################################
    # Step 3: Reconstruct the bands based on the interpolation results
    ###############################################################################
    eband, vvband, cband = my_boltztrap2.fite.getBTPbands(
        equivalences, coeffs, lattvec
    )

    ###############################################################################
    # Step 4: Compute the density of states
    ###############################################################################

    dose, dos, vvdos, cdos = my_boltztrap2.bandlib.BTPDOS(
        eband, vvband, erange=(OFFSET - 0.1, np.max(eband)), npts=2000
    )

    volume = np.linalg.det(data.get_lattvec())

    ###############################################################################
    # Step 5: Compute the electron count, the conductivity and the Seebeck
    # coefficient. Take a value of 1e-14 s for the
    # uniform electronic lifetime.
    ###############################################################################
    tau = 1e-14
    # Define the temperatures and chemical potentials we are interested in
    Tr = np.array([500.0])
    margin = 9.0 * BOLTZMANN * Tr.max()
    mur_indices = np.logical_and(
        dose > dose.min() + margin, dose < dose.max() - margin
    )
    mur = dose[mur_indices]

    # Obtain the Fermi integrals required to get the Onsager coefficients
    N, L0, L1, L2, Lm11 = my_boltztrap2.bandlib.fermiintegrals(
        dose, dos, vvdos, mur=mur, Tr=Tr
    )
    # Translate those into Onsager coefficients
    sigma, seebeck, kappa, Hall = my_boltztrap2.bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, volume
    )
    # Rescale the carrier count into a volumetric density in cm**(-3)
    N = -N[0, ...] / (volume / (Meter / 100.0) ** 3)
    # Obtain the scalar conductivity and Seebeck coefficient
    sigma = tau * sigma[0, ...].trace(axis1=1, axis2=2) / 3.0
    seebeck = seebeck[0, ...].trace(axis1=1, axis2=2) / 3.0
    # Compute the scalar power factor
    P = sigma * seebeck * seebeck
    # Transform these quantities to more convenient units
    sigma *= 1e-3  # kS / m
    seebeck *= 1e6  # microvolt / K
    P *= 1e6  # microwatt / m / K**2

    # Create and show the plot
    plt.figure(figsize=(6, 4))
    ax = plt.gca()
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    plt.plot(N, sigma, label=r"$\sigma\;\left(\mathrm{kS\,m^{-1}}\right)$")
    plt.plot(
        N,
        np.abs(seebeck),
        label=r"$\left\vert S\right\vert\;\left(\mathrm{\mu V\,K^{-1}}\right)$",
    )
    plt.plot(
        N,
        P,
        label=r"$S^2\,\sigma\;\left(\mathrm{\mu W\,m^{-1}\,K^{-2}}\right)$",
    )

    # Do the same thing for the analytic DOS
    analytic_dos = (
        (2 * EFFM) ** 1.5 * np.sqrt(dose - OFFSET) / (4 * np.pi**2)
    )
    analytic_dos[np.isnan(analytic_dos)] = 0.0

    analytic_vvdos = np.zeros(np.shape(vvdos))
    analytic_vvdos[0, 0] = (
        (2 * EFFM) ** 0.5 * (dose - OFFSET) ** 1.5 / (3 * np.pi**2)
    )
    analytic_vvdos[1, 1] = (
        (2 * EFFM) ** 0.5 * (dose - OFFSET) ** 1.5 / (3 * np.pi**2)
    )
    analytic_vvdos[2, 2] = (
        (2 * EFFM) ** 0.5 * (dose - OFFSET) ** 1.5 / (3 * np.pi**2)
    )
    analytic_vvdos[np.isnan(analytic_vvdos)] = 0.0

    N, L0, L1, L2, Lm11 = my_boltztrap2.bandlib.fermiintegrals(
        dose, analytic_dos, analytic_vvdos, mur=mur, Tr=Tr
    )
    sigma, seebeck, kappa, Hall = my_boltztrap2.bandlib.calc_Onsager_coefficients(
        L0, L1, L2, mur, Tr, 1.0
    )
    N = -N[0, ...] / (1.0 / (Meter / 100.0) ** 3)
    sigma = tau * sigma[0, ...].trace(axis1=1, axis2=2) / 3.0
    seebeck = seebeck[0, ...].trace(axis1=1, axis2=2) / 3.0
    P = sigma * seebeck * seebeck
    sigma *= 1e-3  # kS / m
    seebeck *= 1e6  # microvolt / K
    P *= 1e6  # microwatt / m / K**2
    plt.plot(N, sigma, "C0:")
    plt.plot(N, np.abs(seebeck), "C1:")
    plt.plot(N, P, "C2:")

    ax.set_xscale("log")
    plt.xlim(xmin=1e18, xmax=1e21)
    plt.ylim(ymin=0.0, ymax=1000)
    plt.xlabel(r"$N\;\left(\mathrm{cm^{-3}}\right)$", fontsize=16)
    plt.legend(loc="best", fontsize=12)
    plt.tight_layout()
    plt.savefig("IMG_3DHEG2.pdf")
    plt.show()
