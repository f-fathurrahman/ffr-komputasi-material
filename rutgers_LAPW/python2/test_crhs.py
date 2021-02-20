import numpy as np
from solve_scheq import CRHS

RmaxAtom = 10.0  # The end of the radial mesh (maximum r)
Nratom = 3001    # Number of points in radial mesh
R0 = np.linspace(1e-10, RmaxAtom, Nratom) 
Veff = -np.ones(len(R0), dtype=float)/R0

E = -1.0
l = 1
rhs = CRHS(E, l, R0, Veff)

for ir in range(4):
    print("%5d %18.10f %18.10e %18.10e" % (ir+1, R0[ir], Veff[ir], rhs[ir]))

ir = Nratom-2
print("%5d %18.10f %18.10e %18.10e" % (ir+1, R0[ir], Veff[ir], rhs[ir]))

ir = Nratom-1
print("%5d %18.10f %18.10e %18.10e" % (ir+1, R0[ir], Veff[ir], rhs[ir]))
