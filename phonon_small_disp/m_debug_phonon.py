import numpy as np
from math import pi, sqrt
import numpy.linalg as la
import numpy.fft as fft
import ase.units as units

def debug_phonon_read(phonon, method="Frederiksen", symmetrize=3, acoustic=True,
            cutoff=None, born=False, **kwargs):

    print("Enter debug_phonon_read")

    method = method.lower()
    assert method in ["standard", "frederiksen"]
    print("method = ", method)

    if cutoff is not None:
        cutoff = float(cutoff)

    # Read Born effective charges and optical dielectric tensor
    if born:
        phonon.read_born_charges(**kwargs)

    # Number of atoms
    natoms = len(phonon.indices)
    print("natoms = ", natoms)

    # Number of unit cells
    N = np.prod(phonon.supercell)
    print("Number of unit cells = ", N)

    # Matrix of force constants as a function of unit cell index in units
    # of eV / Ang**2
    C_xNav = np.empty((natoms * 3, N, natoms, 3), dtype=float)

    # Loop over all atomic displacements and calculate force constants
    for i, a in enumerate(phonon.indices):
        for j, v in enumerate("xyz"):
            # Atomic forces for a displacement of atom a in direction v
            # basename = "%s.%d%s" % (phonon.name, a, v)
            basename = "%d%s" % (a, v)
            fminus_av = phonon.cache[basename + "-"]["forces"]
            fplus_av = phonon.cache[basename + "+"]["forces"]

            if method == "frederiksen":
                fminus_av[a] -= fminus_av.sum(0)
                fplus_av[a] -= fplus_av.sum(0)

            # Finite difference derivative
            C_av = fminus_av - fplus_av
            C_av /= 2 * phonon.delta

            # Slice out included atoms
            C_Nav = C_av.reshape((N, len(phonon.atoms), 3))[:, phonon.indices]
            index = 3 * i + j
            C_xNav[index] = C_Nav

    # Make unitcell index the first and reshape
    C_N = C_xNav.swapaxes(0, 1).reshape((N,) + (3 * natoms, 3 * natoms))

    # Cut off before symmetry and acoustic sum rule are imposed
    if cutoff is not None:
        phonon.apply_cutoff(C_N, cutoff)

    # Symmetrize force constants
    if symmetrize:
        for i in range(symmetrize):
            # Symmetrize
            C_N = phonon.symmetrize(C_N)
            # Restore acoustic sum-rule
            if acoustic:
                phonon.acoustic(C_N)
            else:
                break

    # Store force constants and dynamical matrix
    phonon.C_N = C_N
    phonon.D_N = C_N.copy()

    # Add mass prefactor
    m_a = phonon.atoms.get_masses()
    phonon.m_inv_x = np.repeat(m_a[phonon.indices]**-0.5, 3)
    M_inv = np.outer(phonon.m_inv_x, phonon.m_inv_x)
    for D in phonon.D_N:
        D *= M_inv



#
def get_band_structure(phonon, path, modes=False, born=False, verbose=True):
    omega_kl = band_structure(phonon, path.kpts, modes, born, verbose)
    if modes:
        assert 0
        omega_kl, modes = omega_kl

    from ase.spectrum.band_structure import BandStructure
    bs = BandStructure(path, energies=omega_kl[None])
    return bs


def band_structure(phonon, path_kc, modes=False, born=False, verbose=True):
    """Calculate phonon dispersion along a path in the Brillouin zone.

    The dynamical matrix at arbitrary q-vectors is obtained by Fourier
    transforming the real-space force constants. In case of negative
    eigenvalues (squared frequency), the corresponding negative frequency
    is returned.

    Frequencies and modes are in units of eV and Ang/sqrt(amu),
    respectively.

    Parameters:

    path_kc: ndarray
        List of k-point coordinates (in units of the reciprocal lattice
        vectors) specifying the path in the Brillouin zone for which the
        dynamical matrix will be calculated.
    modes: bool
        Returns both frequencies and modes when True.
    born: bool
        Include non-analytic part given by the Born effective charges and
        the static part of the high-frequency dielectric tensor. This
        contribution to the force constant accounts for the splitting
        between the LO and TO branches for q -> 0.
    verbose: bool
        Print warnings when imaginary frequncies are detected.

    """

    assert phonon.D_N is not None
    if born:
        assert phonon.Z_avv is not None
        assert phonon.eps_vv is not None

    # Dynamical matrix in real-space
    D_N = phonon.D_N

    # Lists for frequencies and modes along path
    omega_kl = []
    u_kl = []

    # Reciprocal basis vectors for use in non-analytic contribution
    reci_vc = 2 * pi * la.inv(phonon.atoms.cell)
    # Unit cell volume in Bohr^3
    vol = abs(la.det(phonon.atoms.cell)) / units.Bohr**3

    for q_c in path_kc:

        # Add non-analytic part
        if born:
            # q-vector in cartesian coordinates
            q_v = np.dot(reci_vc, q_c)
            # Non-analytic contribution to force constants in atomic units
            qdotZ_av = np.dot(q_v, phonon.Z_avv).ravel()
            C_na = (4 * pi * np.outer(qdotZ_av, qdotZ_av) /
                    np.dot(q_v, np.dot(phonon.eps_vv, q_v)) / vol)
            phonon.C_na = C_na / units.Bohr**2 * units.Hartree
            # Add mass prefactor and convert to eV / (Ang^2 * amu)
            M_inv = np.outer(phonon.m_inv_x, phonon.m_inv_x)
            D_na = C_na * M_inv / units.Bohr**2 * units.Hartree
            phonon.D_na = D_na
            D_N = phonon.D_N + D_na / np.prod(phonon.supercell)

        # if np.prod(phonon.N_c) == 1:
        #
        #     q_av = np.tile(q_v, len(phonon.indices))
        #     q_xx = np.vstack([q_av]*len(phonon.indices)*3)
        #     D_m += q_xx

        # Evaluate fourier sum
        D_q = compute_dynamical_matrix(phonon, q_c, D_N)

        if modes:
            omega2_l, u_xl = la.eigh(D_q, UPLO='U')
            # Sort eigenmodes according to eigenvalues (see below) and
            # multiply with mass prefactor
            u_lx = (phonon.m_inv_x[:, np.newaxis] *
                    u_xl[:, omega2_l.argsort()]).T.copy()
            u_kl.append(u_lx.reshape((-1, len(phonon.indices), 3)))
        else:
            omega2_l = la.eigvalsh(D_q, UPLO='U')

        # Sort eigenvalues in increasing order
        omega2_l.sort()
        # Use dtype=complex to handle negative eigenvalues
        omega_l = np.sqrt(omega2_l.astype(complex))

        # Take care of imaginary frequencies
        if not np.all(omega2_l >= 0.):
            indices = np.where(omega2_l < 0)[0]

            if verbose:
                print('WARNING, %i imaginary frequencies at '
                        'q = (% 5.2f, % 5.2f, % 5.2f) ; (omega_q =% 5.3e*i)'
                        % (len(indices), q_c[0], q_c[1], q_c[2],
                            omega_l[indices][0].imag))

            omega_l[indices] = -1 * np.sqrt(np.abs(omega2_l[indices].real))

        omega_kl.append(omega_l.real)

    # Conversion factor: sqrt(eV / Ang^2 / amu) -> eV
    s = units._hbar * 1e10 / sqrt(units._e * units._amu)
    omega_kl = s * np.asarray(omega_kl)

    if modes:
        return omega_kl, np.asarray(u_kl)

    return omega_kl


def compute_dynamical_matrix(phonon, q_scaled: np.ndarray, D_N: np.ndarray):
    """ Computation of the dynamical matrix in momentum space D_ab(q).
        This is a Fourier transform from real-space dynamical matrix D_N
        for a given momentum vector q.

    q_scaled: q vector in scaled coordinates.

    D_N: the dynamical matrix in real-space. It is necessary, at least
            currently, to provide this matrix explicitly (rather than use
            self.D_N) because this matrix is modified by the Born charges
            contributions and these modifications are momentum (q) dependent.

    Result:
        D(q): two-dimensional, complex-valued array of
                shape=(3 * natoms, 3 * natoms).
    """
    # Evaluate fourier sum
    R_cN = phonon._lattice_vectors_array
    phase_N = np.exp(-2.j * pi * np.dot(q_scaled, R_cN))
    D_q = np.sum(phase_N[:, np.newaxis, np.newaxis] * D_N, axis=0)
    return D_q