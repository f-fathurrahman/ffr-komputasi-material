from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi
import scipy

from atom_cmpb import *
from Numerov import Numerov
from solve_scheq import CRHS

from find_core_states import find_core_states

# Some parameters
RmaxAtom = 10.0  # The end of the radial mesh (maximum r)
Nratom = 3001    # Number of points in radial mesh

Z = 29
# Core specification
core = [3,2,0]

print("Z           = %d" % Z)
print("core config = ", core)

# Radial mesh, equally spaced
R0 = np.linspace(1e-10, RmaxAtom, Nratom) 
for ir in range(5):
    print("%d %18.10f" % (ir, R0[ir]))

# Inverse radial mesh (inverted in sequence)
Ra = R0[::-1]
for ir in range(5):
    print("%d %18.10f" % (ir, Ra[ir]))

# We add one more state to core to get atomic states
catm = [c + 1 for c in core]

Veff = -Z*np.ones(len(Ra), dtype=float)/Ra

(coreRho, coreE, coreZ, states) = find_core_states(catm, Ra, Veff, Z, fraction=4.0)

print("coreE = ")
print(coreE)

l_label = ["s", "p", "d", "f"]

print("states = ")
for ist in range(len(states)):
    n, l, E = states[ist]
    print("%d%s = %18.10f" % (n+1, l_label[l], E ))

# Sorts them according to energy
states.sort(atom_cmpb)
print()
print("After sort:")
for ist in range(len(states)):
    n, l, E = states[ist]
    print("%d%s = %18.10f" % (n+1, l_label[l], E ))
