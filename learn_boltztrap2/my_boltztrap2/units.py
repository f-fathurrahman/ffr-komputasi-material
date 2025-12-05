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

import math

# Values retrieved from physics.nist.gov/constants on 2017-08-04
BOLTZMANN_SI = 1.38064852e-23
FINE_STRUCT = 7.2973525664e-3
AVOGADRO = 6.022140857e23  # [1/mol]
Clight_SI = 299792458.0  # [m/s]
a0_SI = 5.2917721067e-11  # [m]
me_SI = 9.10938356e-31  # [kg]
qe_SI = 1.6021766208e-19  # Electron Charge [C]

Clight = 1.0 / FINE_STRUCT
EPSILON0 = 1.0 / (4.0 * math.pi)
MU0 = 1.0 / (EPSILON0 * Clight**2)
MUB = 1.0 / 2.0

Meter = 1.0 / a0_SI
Kilogram = 1.0 / me_SI
Coulomb = 1.0 / qe_SI
Second = Clight_SI / Clight * Meter

Newton = Kilogram * Meter / Second**2
Joule = Newton * Meter
Volt = Joule / Coulomb
Ampere = Coulomb / Second
Ohm = Volt / Ampere
Siemens = Ampere / Volt

hbar_SI = 1 / (Joule * Second)
Angstrom = Meter * 1e-10
eV = Volt
BOLTZMANN = BOLTZMANN_SI * Joule
Ha = 1.0
