# -*- coding: utf-8 -*
#    my_boltztrap2, a program for interpolating band structures and calculating
#    semi-classical transport coefficients.
#    Copyright (C) 2017-2025 Georg K. H. Madsen <georg.madsen@tuwien.ac.at>
#    Copyright (C) 2017-2025 Jesús Carrete <jesus.carrete.montana@tuwien.ac.at>
#    Copyright (C) 2017-2025 Matthieu J. Verstraete <matthieu.verstraete@ulg.ac.be>
#    Copyright (C) 2018-2019 Genadi Naydenov <gan503@york.ac.uk>
#    Copyright (C) 2020 Gavin Woolman <gwoolma2@staffmail.ed.ac.uk>
#    Copyright (C) 2020 Roman Kempt <roman.kempt@tu-dresden.de>
#    Copyright (C) 2022 Robert Stanton <stantor@clarkson.edu>
#    Copyright (C) 2024 Haoyu (Daniel) Yang <yanghaoyu97@outlook.com>
#
#    This file is part of my_boltztrap2.
#
#    my_boltztrap2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    my_boltztrap2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with my_boltztrap2. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy as sp
import scipy.integrate
import scipy.optimize
import scipy.signal

from my_boltztrap2.fd import *
from my_boltztrap2.fd import _FD_XMAX, _FD_XMAX_GAP
from my_boltztrap2.units import *


def _suggest_nbins(data, erange, Tmin=None):
    """Suggest a number of bins for a histogram based on numpy's criterion.

    The number returned is the maximum of the prescriptions of the Sturges
    and Freedman-Draconis rules, which should be roughly equivalent to what
    numpy.histogram would use when nbins="auto", or the minimum number of
    binds such that the energy difference between adjacent bins is at least
    0.5 * kB * T at the lowest temperature studied (whatever is higher).
    If Tmin is None, the temperature-based criterion is not used.
    """
    # Reimplementation of the two parts of the NumPy "auto" criterion.
    flat = data.ravel()
    sturges = np.log2(flat.size) + 1.0
    iqr = np.percentile(flat, 75.0) - np.percentile(flat, 25.0)
    h = 2.0 * iqr / flat.size ** (1.0 / 3.0)
    fd = (erange[1] - erange[0]) / h
    if Tmin is None:
        return int(np.ceil(max(sturges, fd)))
    # Temperature-based criterion.
    bytemp = 2.0 * (erange[1] - erange[0]) / (BOLTZMANN * Tmin)
    return int(np.ceil(max(sturges, fd, bytemp)))


def DOS(eigs, erange=None, npts=None, weights=None, Tmin=None):
    """Compute the density of states.

    Args:
        eband: (nbands, nkpoints) array with the band energies
        erange: 2-tuple with the minimum and maximum energies to be considered.
            If its value is None, take the minimum and maximum band energies.
        npts: number of bins to include in the histogram. If omitted,
            _suggest_nbins will be called to obtain an estimate.
        weights: array with the same shape as eband to be used as the weights.
        Tmin: minimum temperature explored in the calculation. If provided, and
            if npts is not none, it will be passed to _suggest_nbins.

    Returns:
        Two 1D numpy arrays of the same size with the bin energies and the DOS,
        respectively.
    """
    nkpt, nband = np.shape(eigs)
    if erange is None:
        erange = (eigs.min(), eigs.max())
    if npts is None:
        npts = _suggest_nbins(eigs, erange, Tmin)
    pip = np.histogram(eigs, npts, weights=weights, range=erange)
    npts = pip[1].size - 1
    tdos = np.zeros((2, npts), dtype=float)
    tdos[1] = pip[0] / float(nkpt) / ((erange[1] - erange[0]) / npts)
    tdos[0] = 0.5 * (pip[1][:-1] + pip[1][1:])
    return tdos


def lambda_to_tau(vvband, lambdas):
    """Translate a set of mean free paths into a set of lifetimes.

    Args:
        vvband: 4D numpy array containing the outer product of the group
            velocity of each electron mode with itself. It must have shape
            (nbands, 3, 3, nkpoints).
        lambdas: 2D numpy array with a mean free path for each electron mode.
            It must have shape (nbands, nkpoints).

    Returns:
        A vector with the same shape as the lambdas argument, containing the
        lifetime for each electron mode.
    """
    if not isinstance(lambdas, np.ndarray) or lambdas.ndim != 2:
        raise ValueError("lambdas must be a 2D numpy array")
    if (
        not isinstance(vvband, np.ndarray)
        or vvband.ndim != 4
        or vvband.shape != (lambdas.shape[0], 3, 3, lambdas.shape[1])
    ):
        raise ValueError("vvband has an invalid shape")
    # Extract the norm of each group velocity
    scalarv2 = np.trace(vvband, axis1=1, axis2=2)
    scalarv = np.sqrt(scalarv2)
    # Perform the transformation while avoiding divisions by zero.
    nruter = np.zeros_like(lambdas)
    nonzero = scalarv2 >= 1e-8
    nruter[nonzero] = lambdas[nonzero] / scalarv[nonzero]
    return nruter


def BTPDOS(
    eband,
    vvband,
    cband=None,
    erange=None,
    npts=None,
    scattering_model="uniform_tau",
    Tmin=None,
):
    """Compute the DOS, transport DOS and "curvature DOS".

    The transport DOS is weighted by the outer product of the group velocity
    with itself, and by the relaxation time. The "curvature DOS" is weighted by
    the curvature and the relaxation time.

    In order to use a custom mean free path instead of a custom lifetime, use
    the auxiliary function lambda_to_tau().

    Args:
        eband: (nbands, nkpoints) array with the band energies
        vvband: (nbands, 3, 3, nkpoints) array with the outer product of each
            group velocity with itself.
        cband: None or (nbands, 3, 3, 3, nkpoints) array with the curvatures
        erange: range of energies for the DOS. If not provided, it will be
            automatically determined by the DOS function.
        npts: number of bins in the histogram used to determine the DOS. If not
            provided, the DOS function will take a conservative guess.
        scattering model: model to be used for the electron lifetimes. The
            following choices are available:
                - "uniform_tau": uniform lifetime for all carriers
                - "uniform_lambda": uniform mean free path for all carriers
                - A 2d array with the same shape as eband, with a scattering
                  rate for each electron mode.
        Tmin: minimum temperature explored in the calculation. If provided, it
            will be passed to DOS.

    Returns:
        Four arrays containing the bin energies, the DOS, the transport dos and
        the "curvature DOS". If cband is none, the last element of the returned
        value will also be none. The sizes of the returned arrays are (npts, ),
        (npts,), (3, 3, npts) and (3, 3, 3, npts).
    """
    dos = DOS(eband.T, erange=erange, npts=npts, Tmin=Tmin)
    npts = dos[0].size
    iu0 = np.array(np.triu_indices(3)).T
    vvdos = np.zeros((3, 3, npts))
    multpl = np.ones_like(eband)
    if isinstance(scattering_model, str) and scattering_model == "uniform_tau":
        pass
    elif (
        isinstance(scattering_model, str)
        and scattering_model == "uniform_lambda"
    ):
        multpl = lambda_to_tau(vvband, multpl)
    elif isinstance(scattering_model, np.ndarray):
        if scattering_model.shape != eband.shape:
            raise ValueError(
                "scattering_model and ebands must have the same shape"
            )
        multpl = scattering_model
    else:
        raise ValueError("unknown scattering model")
    for i, j in iu0:
        weights = vvband[:, i, j, :] * multpl
        vvdos[i, j] = DOS(
            eband.T, weights=weights.T, erange=erange, npts=npts
        )[1]
    il1 = np.tril_indices(3, -1)
    iu1 = np.triu_indices(3, 1)
    vvdos[il1[0], il1[1]] = vvdos[iu1[0], iu1[1]]

    if cband is None:
        cdos = None
    else:
        cdos = np.zeros((3, 3, 3, npts))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    weights = cband[:, i, j, k, :] * multpl * multpl
                    cdos[i, j, k] = DOS(
                        eband.T, weights=weights.T, erange=erange, npts=npts
                    )[1]
    return dos[0], dos[1], vvdos, cdos


def smoothen_DOS(e, dos, Tsmooth):
    """Return a smoothed version of the DOS.

    The smoothed DOS is obtained by convolving the original with the derivative
    of the FD distribution at a given temperature. Note that boundary effects
    are not corrected for.

    Args:
        e: array of energies
        dos: original ("rough") density of states
        Tsmooth: temperature parameter for the smoothing Fermi-Dirac kernel

    Returns:
        The smoothed DOS, with the same shape as the "rough" one.
    """
    de = e[1] - e[0]
    kernel = np.abs(dFDde(e - e.mean(), 0.0, BOLTZMANN * Tsmooth))
    kernel = np.trim_zeros(kernel)
    if len(kernel) == 0:
        raise ValueError(
            "no non-zero values in DOS convolution kernel. Grid is probably too coarse"
        )
    nruter = sp.signal.convolve(dos, kernel, mode="same")
    nruter *= dos.sum() / nruter.sum()
    return nruter


def smoothen_DOS_direct(e, dos, Tsmooth, ep):
    """Return a smoothed version of the DOS evaluated at the energies in ep.

    The smoothed DOS is obtained by convolving the original with the derivative
    of the FD distribution at a given temperature. The convolution is evaluated
    directly.

    Args:
        e: array of input energies
        dos: original ("rough") density of states
        Tsmooth: temperature parameter for the smoothing Fermi-Dirac kernel
        ep: array of energies at which to evaluate the smooth dos

    Returns:
        The smoothed DOS, evaluated at the energies in ep.
    """
    global _FD_XMAX

    try:
        # Only available in NumPy >= 2.0.
        integrator = np.trapezoid
    except AttributeError:
        integrator = np.trapz
    integrator = np.trapz

    norm = sp.integrate.quad(
        lambda x: dFDde(x, 0.0, BOLTZMANN * Tsmooth),
        -_FD_XMAX * BOLTZMANN * Tsmooth,
        _FD_XMAX * BOLTZMANN * Tsmooth,
    )[0]
    nruter = np.zeros_like(ep)
    for i, mu in np.ndenumerate(ep):
        d = dFDde(mu - e, 0.0, BOLTZMANN * Tsmooth)
        nruter[i] = -integrator(d * dos, e)
    nruter /= norm
    return nruter


def calc_N(epsilon, dos, mu, T, dosweight=2.0):
    """Compute the electron count by integrating over the DOS.

    Args:
        epsilon: array of energies at which the DOS is available
        dos: density of states
        mu: single value of the chemical potential
        T: single value of the temperature
        dosweight: maximum occupancy of an electron mode

    Returns:
        A single value of the electron count for the given values of mu and T.
    """
    if T == 0.0:
        occ = np.where(epsilon < mu, 1.0, 0.0)
        occ[epsilon == mu] = 0.5
    else:
        kBT = T * BOLTZMANN
        occ = FD(epsilon, mu, kBT)
    de = epsilon[1] - epsilon[0]
    return -dosweight * (dos * occ).sum() * de


def solve_for_mu(
    epsilon, dos, N0, T, dosweight=2.0, refine=False, try_center=False
):
    """Estimate the chemical potential required to have N0 electrons.

    First, a value of mu is obtained that yields a number of electrons as
    close as possible to N0. If mu falls in a wide gap (relative to kB * T),
    then mu0 is moved to the center of the gap.

    Args:
        epsilon: array of energies at which the DOS is available
        dos: density of states
        N0: number of valence electrons in the compound
        T: single value of the temperature
        dosweight: maximum occupancy of an electron mode
        refine: if True, do not restrict mu to the values in epsilon
        try_center: if True, snap the value of mu to the center of large gaps

    Returns:
        An estimate of the intrinsic chemical potential at a given temperature.
    """
    global _FD_XMAX_GAP
    # Compute the distance to electroneutrality at each sampled energy
    delta = np.empty_like(epsilon)
    for i, e in enumerate(epsilon):
        delta[i] = calc_N(epsilon, dos, e, T, dosweight) + N0
    delta = np.abs(delta)
    # Find the position optimizing this distance
    pos = np.abs(delta).argmin()
    mu = epsilon[pos]
    center = False
    # Check if mu0 falls in a gap
    if dos[pos] == 0.0:
        # If it does, search for the edges of the gap
        lpos = -1
        hpos = -1
        for i in range(pos, -1, -1):
            if dos[i] != 0.0:
                lpos = i
                break
        for i in range(pos, dos.size):
            if dos[i] != 0.0:
                hpos = i
                break
        if -1 in (lpos, hpos):
            raise ValueError("mu0 lies outside the range of band energies")
        # If mu0 is in a gap, and far enough from the edges, move it to the
        # center of the gap if requested.
        hepsilon = epsilon[hpos]
        lepsilon = epsilon[lpos]
        kBT = T * BOLTZMANN
        if (
            try_center
            and min(hepsilon - mu, mu - lepsilon) >= _FD_XMAX_GAP * kBT / 2.0
        ):
            pos = int(round(0.5 * (lpos + hpos)))
            mu = epsilon[pos]
            center = True
    if refine:
        # Either move to the exact center of the gap, or look for the best-fit
        # value of mu in the selected histogram bin.
        if center:
            mu = 0.5 * (lepsilon + hepsilon)
        else:
            residual = calc_N(epsilon, dos, mu, T, dosweight) + N0
            if np.isclose(residual, 0):
                lpos = pos
                hpos = pos
            elif residual > 0:
                lpos = pos
                hpos = min(pos + 1, epsilon.size - 1)
            else:
                lpos = max(0, pos - 1)
                hpos = pos
            if hpos != lpos:
                lmu = epsilon[lpos]
                hmu = epsilon[hpos]

                def calc_abs_residual(muarg):
                    """Auxiliary function to be minimized."""
                    return abs(calc_N(epsilon, dos, muarg, T, dosweight) + N0)

                result = sp.optimize.minimize_scalar(
                    calc_abs_residual, bounds=(lmu, hmu), method="bounded"
                )
                mu = result.x
    return mu


def calc_cv(epsilon, dos, mur, Tr, dosweight=2.0):
    """Compute the electronic contribution to the heat capacity.

    Args:
        epsilon: array of energies at which the DOS is available
        dos: density of states
        N0: number of valence electrons in the compound
        mur: array of chemical potential values
        Tr: array of temperature values
        dosweight: maximum occupancy of an electron mode

    Returns:
        An array with the same shape as the outer product of Tr and mur,
        containing the electronic contribution to the heat capacity for each
        temperature and chemical potential, in SI units.
    """
    kBTr = np.array(Tr) * BOLTZMANN
    nT = len(Tr)
    nmu = len(mur)
    nruter = np.empty((nT, nmu))
    de = epsilon[1] - epsilon[0]
    for iT, kBT in enumerate(kBTr):
        for imu, mu in enumerate(mur):
            nruter[iT, imu] = (
                -(
                    dosweight
                    * dos
                    * dFDde(epsilon, mu, kBT)
                    * (epsilon - mu) ** 2
                    / kBT
                ).sum()
                * de
            )
    return BOLTZMANN_SI * nruter


def fermiintegrals(epsilon, dos, sigma, mur, Tr, dosweight=2.0, cdos=None):
    """Compute the moments of the FD distribution over the band structure.

    Args:
        epsilon: array of energies at which the DOS is available
        dos: density of states
        sigma: transport DOS
        mur: array of chemical potential values
        Tr: array of temperature values
        dosweight: maximum occupancy of an electron mode
        cdos: "curvature DOS" if available

    Returns:
        Five numpy arrays, namely:
        1. An (nT, nmu) array with the electron counts for each temperature and
           each chemical potential.
        2. An (nT, nmu, 3, 3) with the integrals of the 3 x 3 transport DOS
           over the band structure taking the occupancies into account.
        3. An (nT, nmu, 3, 3) with the first moment of the 3 x 3 transport DOS
           over the band structure taking the occupancies into account.
        4. An (nT, nmu, 3, 3) with the second moment of the 3 x 3 transport DOS
           over the band structure taking the occupancies into account.
        5. If the cdos argument is provided, an (nT, nmu, 3, 3, 3) with the
           integrals of the 3 x 3 x 3 "curvature DOS" over the band structure
           taking the occupancies into account.
        where nT and nmu are the sizes of Tr and mur, respectively.
    """
    kBTr = np.array(Tr) * BOLTZMANN
    iu0 = np.triu_indices(3)
    ido = np.tril_indices(3, -1)
    iu1 = np.triu_indices(3, 1)
    nT = len(Tr)
    nmu = len(mur)
    N = np.empty((nT, nmu))
    L0 = np.empty((nT, nmu, 3, 3))
    L1 = np.empty((nT, nmu, 3, 3))
    L2 = np.empty((nT, nmu, 3, 3))
    if cdos is not None:
        L11 = np.empty((nT, nmu, 3, 3, 3))
    else:
        L11 = None
    de = epsilon[1] - epsilon[0]
    for iT, kBT in enumerate(kBTr):
        for imu, mu in enumerate(mur):
            N[iT, imu] = -(dosweight * dos * FD(epsilon, mu, kBT)).sum() * de
            int0 = -dosweight * dFDde(epsilon, mu, kBT)
            intn = int0 * sigma
            L0[iT, imu] = intn.sum(axis=2) * de
            intn *= epsilon - mu
            L1[iT, imu] = -intn.sum(axis=2) * de
            intn *= epsilon - mu
            L2[iT, imu] = intn.sum(axis=2) * de
            if cdos is not None:
                cint = int0 * cdos
                L11[iT, imu] = -cint.sum(axis=3) * de
    return N, L0, L1, L2, L11


def calc_Onsager_coefficients(L0, L1, L2, mur, Tr, vuc, Lm11=None):
    """Compute a set of Onsager coefficients based on Fermi-Dirac integrals.

    The input parameters will most likely come from fermiintegrals().

    Args:
        L0: an array with the integrals of the 3 x 3 transport DOS over the
            band structure taking the occupancies into account, for each value
            of the temperature and the chemical  potential.
        L1: an array with the first moment of the 3 x 3 transport DOS over the
            band structure taking the occupancies into account, for each value
            of the temperature and the chemical  potential.
        L2: an array with the second moment of the 3 x 3 transport DOS over the
            band structure taking the occupancies into account, for each value
            of the temperature and the chemical  potential.
        mur: array of chemical potential values
        Tr: array of temperature values
        vuc: volume of the unit cell
        Lm11: optional array with the integrals of the 3 x 3 x 3 "curvature
            DOS" over the band structure taking the occupancies into account,
            for each value of the temperature and the chemical  potential.

    Returns:
        A set of four numpy arrays:
        1. A (nT, nmu, 3, 3) array with the 3 x 3 conductivity tensor for each
           combination of temperature and chemical potential.
        2. A (nT, nmu, 3, 3) array with the 3 x 3 Seebeck coefficient tensor
           for each combination of temperature and chemical potential.
        3. A (nT, nmu, 3, 3) array with the 3 x 3 charge carrier contribution
           to the thermal conductivity at zero current, for each combination of
           temperature and chemical potential.
        4. If the Lm11 argument is provided, a (nT, nmu, 3, 3, 3) array with
           the 3 x 3 x 3 Hall tensor, for each combination of temperature and
           chemical potential.
        where nT and nmu are the sizes of Tr and mur, respectively.
    """
    L11 = L0 / (Siemens / (Meter * Second)) / vuc
    L12 = (
        L1
        / Tr[:, None, None, None]
        / (Volt * Siemens / (Meter * Second))
        / vuc
    )
    L22 = (
        L2
        / Tr[:, None, None, None]
        / (Volt * Joule * Siemens / (Meter * Second * Coulomb))
        / vuc
    )
    if Lm11 is not None:
        L111 = Lm11 / (Siemens**2 * Meter / (Second**3 * Ampere)) / vuc
    seebeck = np.empty_like(L11)
    kappa = np.empty_like(L22)
    if Lm11 is not None:
        Hall = np.empty_like(L111)
    else:
        Hall = None
    for iT in range(L11.shape[0]):
        for imu in range(L11.shape[1]):
            pL11 = np.linalg.pinv(L11[iT, imu])
            seebeck[iT, imu, :, :] = pL11 @ L12[iT, imu, :, :]
            kappa[iT, imu, :, :] = (
                L22[iT, imu, :, :]
                - Tr[iT]
                * L11[iT, imu, :, :]
                @ seebeck[iT, imu, :, :]
                @ seebeck[iT, imu, :, :]
            )
            if Lm11 is not None:
                Hall_ = np.zeros((3, 3, 3))
                for k in range(3):
                    Hall_[:, :, k] = pL11 @ L111[iT, imu, :, :, k] @ pL11
                Hall[iT, imu] = Hall_
    return (L11, seebeck, kappa, Hall)


def apply_scissor(epsilon, dos, N0, eband, desired_gap, dosweight=2.0):
    """Shift all the conduction bands to achieve the desired value of the gap.

    Do not forget to readjust the Fermi level after shifting the bands!

    Args:
        epsilon: array of energies at which the DOS is available
        dos: density of states
        N0: number of valence electrons in the compound
        eband: (nbands, nkpoints) array with the band energies
        dosweight: maximum occupancy of an electron mode

    Returns:
        An array with the new band energies.

    Raises:
        ValueError: if there is no gap to begin with or if the Fermi energy
            lies outside of the interval spanned by the bands.
    """
    # Start by finding the current Fermi energy.
    fermi = solve_for_mu(
        epsilon, dos, N0, 0.0, dosweight, refine=False, try_center=False
    )
    # Check if "fermi" lies in a gap, and raise an exception otherwise.
    pos = np.abs(epsilon - fermi).argmin()
    if dos[pos] != 0.0:
        raise ValueError("there is no gap to adjust")
    # Search for the edges of the gap.
    lpos = -1
    hpos = -1
    for i in range(pos, -1, -1):
        if dos[i] != 0.0:
            lpos = i
            break
    for i in range(pos, dos.size):
        if dos[i] != 0.0:
            hpos = i
            break
    if -1 in (lpos, hpos):
        raise ValueError(
            "the Fermi energy lies outside the range of " + "band energies"
        )
    current_gap = epsilon[hpos] - epsilon[lpos]
    delta = desired_gap - current_gap
    # Detect all bands that lie entirely above the gap.
    conduction_mask = eband.min(axis=1) > fermi
    # And shift them. In most cases the goal will be to increase the gap,
    # but in principle it will also work for reducing it.
    nruter = eband.copy()
    nruter[conduction_mask, :] += delta
    return nruter
