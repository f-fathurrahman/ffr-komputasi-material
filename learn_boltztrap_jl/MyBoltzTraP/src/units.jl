# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
# Part of BoltzTraP.jl

# ============================================================================
# Unit Conversion Constants
# ============================================================================
#
# This module provides numeric conversion constants for unit handling.
# Following Python BoltzTraP2's approach, we use explicit numeric constants
# rather than Unitful.jl to ensure exact compatibility with BoltzTraP2 results.
#
# All values are from CODATA 2018.
# ============================================================================

# Length: Bohr ↔ Ångström
const BOHR_TO_ANG = 0.529177210903      # 1 Bohr = 0.529177210903 Å (exact)
const ANG_TO_BOHR = 1 / BOHR_TO_ANG     # 1 Å = 1.8897259886 Bohr

# Energy: Hartree ↔ eV
const HA_TO_EV = 27.211386245988        # 1 Ha = 27.211386245988 eV (CODATA 2018)
const EV_TO_HA = 1 / HA_TO_EV           # 1 eV = 0.03674932217565 Ha

# ============================================================================
# Boltzmann Constant
# ============================================================================

# Boltzmann constant in atomic units (Hartree/Kelvin)
# k_B = 1.380649e-23 J/K (exact, SI 2019)
# 1 Ha = 4.3597447222071e-18 J
# k_B [Ha/K] = 1.380649e-23 / 4.3597447222071e-18 = 3.1668115634438576e-6
const KB_AU = 3.1668115634438576e-6     # Ha/K

# Boltzmann constant in SI units (Joule/Kelvin)
# Using Python BoltzTraP2 value for exact compatibility (NIST 2017)
const KB_SI = 1.38064852e-23            # J/K

# ============================================================================
# Onsager Coefficient Conversion Factors
# ============================================================================
#
# These convert transport coefficients from atomic units to SI units.
# Derivation (using CODATA 2018 constants):
#
#   a0 = 0.529177210903e-10 m (Bohr radius)
#   me = 9.1093837015e-31 kg (electron mass)
#   qe = 1.602176634e-19 C (electron charge, exact)
#   c  = 299792458 m/s (speed of light, exact)
#   α  = 7.2973525693e-3 (fine structure constant)
#
#   c_au = 1/α ≈ 137.036 (speed of light in atomic units)
#   Meter = 1/a0
#   Second = c/c_au * Meter
#   Coulomb = 1/qe
#   Kilogram = 1/me
#   Newton = Kilogram * Meter / Second²
#   Joule = Newton * Meter
#   Volt = Joule / Coulomb
#   Ampere = Coulomb / Second
#   Siemens = Ampere / Volt
#
#   SIGMA_CONV = Siemens / (Meter * Second)           [S/(m·s)]
#   SEEBECK_CONV = Volt * Siemens / (Meter * Second)  [V/K intermediate]
#   KAPPA_CONV = Volt * Joule * Siemens / (Meter * Second * Coulomb)  [W/(m·K·s)]
#
# Pre-computed values (verified against Unitful.jl calculation):

const SIGMA_CONV = 5.2586177987809155e-24    # σ/τ: S/(m·s)
const SEEBECK_CONV = 1.9325063968640632e-25  # Seebeck intermediate: V·S/(m·s)
const KAPPA_CONV = 7.101830018500863e-27     # κ/τ: W/(m·K·s)
