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
import scipy.optimize
import scipy.signal

from my_boltztrap2.units import *

# Threshold for non-zero values of the Fermi-Dirac distribution.
# If |x| >= _FD_XMAX, fFD(x) is approximated by 0.
# _FD_XMAX is defined as the point where fFD(x) == _FD_THRESHOLD
_FD_THRESHOLD = 1e-8

# A much laxer threshold is used to define a wide gap
_FD_THRESHOLD_GAP = 1e-3


def _fd_criterion_gen(threshold):
    """Generate a function to translate the thresholds above into energies.

    Args:
        threshold: threshold for non-zero values of the Fermi-Dirac
            distribution.

    Returns:
        A single-argument function returning fFD(x) - threshold.
    """

    def _fd_criterion(x):
        """Auxiliary function whose root is _FD_XMAX."""
        return 1.0 / (np.exp(x) + 1.0) - threshold

    return _fd_criterion


_FD_XMAX = sp.optimize.newton(_fd_criterion_gen(_FD_THRESHOLD), 0.0)
_FD_XMAX_GAP = sp.optimize.newton(_fd_criterion_gen(_FD_THRESHOLD_GAP), 0.0)

# The same _FD_XMAX is used as the cutoff for evaluationg the derivative of the
# F-D distribution. Note x**n * dFD/dx can take significantly longer to decay
# for large values of n. This is hardly a concern for realistic models of
# scattering and the energy ranges typically afforded by DFT calculations.


def FD(e, mu, kbT):
    """Compute Fermi-Dirac occupancies.

    Args:
        e: array of energies
        mu: single value of the chemical potential
        kBT: thermal energy at the temperature of interest

    Returns:
        An array with the same shape as the "e" argument containing the average
        Fermi-Dirac occupancies for each energy level.
    """
    global _FD_XMAX
    if kbT == 0.0:
        delta = e - mu
        nruter = np.where(delta < 0.0, 1.0, 0.0)
        nruter[np.isclose(delta, 0.0)] = 0.5
    else:
        x = np.asarray((e - mu) / kbT)
        nruter = np.where(x < 0.0, 1.0, 0.0)
        # Note that fFD(0) = 0.5 is always calculated explicitly
        indices = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
        nruter[indices] = 1.0 / (np.exp(x[indices]) + 1.0)
    return nruter


def dFDdx(x):
    """Compute the derivative of the Fermi-Dirac occupancies with respect to
    the recentered and scaled energy (e - mu) / kBT

    Args:
        x: recentered and scaled energies at which to compute the derivative

    Returns:
        An array with the same shape as x containing the derivatives.
    """
    global _FD_XMAX
    x = np.asarray(x)
    nruter = np.zeros_like(x)
    indices = np.logical_and(x > -_FD_XMAX, x < _FD_XMAX)
    c = np.cosh(0.5 * x[indices])
    nruter[indices] = -0.25 / c / c
    return nruter


def dFDde(e, mu, kbT):
    """Compute the derivative of Fermi-Dirac occupancies with respect to
    the energy.

    Args:
        e: array of energies
        mu: single value of the chemical potential
        kBT: thermal energy at the temperature of interest

    Returns:
        An array with the same shape as the "e" argument containing the
        derivatives of the average Fermi-Dirac occupancies for each energy
        level.
    """
    factor = (e - mu) / kbT
    dfde = dFDdx(factor) / kbT
    return dfde
