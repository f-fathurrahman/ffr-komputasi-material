from __future__ import print_function
import numpy as np
from math import exp

from solve_scheq import CRHS
from Numerov import Numerov

RmaxAtom = 10.0  # The end of the radial mesh (maximum r)
Nratom = 3001    # Number of points in radial mesh
R = np.linspace(1e-10, RmaxAtom, Nratom) 
Veff = -np.ones(len(R), dtype=float)/R

E = -1.0
l = 1
rhs = CRHS(E, l, R, Veff)

#for ir in range(4):
#    print("%5d %18.10f %18.10e %18.10e" % (ir+1, R[ir], Veff[ir], rhs[ir]))
#ir = Nratom-2
#print("%5d %18.10f %18.10e %18.10e" % (ir+1, R[ir], Veff[ir], rhs[ir]))
#ir = Nratom-1
#print("%5d %18.10f %18.10e %18.10e" % (ir+1, R[ir], Veff[ir], rhs[ir]))

h = (R[-1]-R[0])/(len(R)-1.)
u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))

for ir in range(4):
    print("%5d %18.10e" % (ir+1, u[ir]))
ir = Nratom-2
print("%5d %18.10e" % (ir+1, u[ir]))
ir = Nratom-1
print("%5d %18.10e" % (ir+1, u[ir]))
