from __future__ import print_function
from Numerov import Numerov
from CRHS import CRHS
import numpy as np
from math import exp

def f_root_sch(Ex, l, R, Veff):
    rhs = CRHS(Ex, l, R, Veff)
    h = (R[-1]-R[0])/(len(R)-1.)
    u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))
    extraplt = u[-2]*(2+h**2*rhs[-2]) - u[-3] # not used?
    return u[-1]

# Some parameters
RmaxAtom = 10.0  # The end of the radial mesh (maximum r)
Nratom = 3001    # Number of points in radial mesh

Z = 1

l = 0
# Radial mesh, equally spaced
R0 = np.linspace(1e-10, RmaxAtom, Nratom) 
for ir in range(5):
    print("%d %18.10f" % (ir, R0[ir]))

# Inverse radial mesh (inverted in sequence)
Ra = R0[::-1]
for ir in range(5):
    print("%d %18.10f" % (ir, Ra[ir]))

Veff = -29*np.ones(len(Ra), dtype=float)/Ra

E = -0.5*Z*Z/(l+1)**2 #- 3.0      # here we starts to look for zero

print("E = ", E)
v0 = f_root_sch(E, l, Ra, Veff)      # starting value
print("v0 = ", v0)

fraction = 4.0
dE = abs(E)/fraction          # the length of the first step
decrease = abs(E)/(abs(E)-dE) # the step will decrease to zero. Check the formula!

#v1 = 0.0
#while E <= 0.0:
#    # and n < core[l]:      # we need ncore[l] bound states
#    #if iterShoot >= NiterShootMax:
#    #    print("shooting method is not converged for n=%d, l=%d" % (n,l))
#    #    print("last E = ", E)
#    #    break
#    E = E + dE
#    print("E = ", E)
#    v1 = f_root_sch(E, l, Ra, Veff)
#    if v1*v0 < 0.0:
#        print("Should find root here")
#        # root solving
#        #Energy = scipy.optimize.brentq(root, E-dE, E, args=(l, R, Veff))
#    dE = dE/decrease
#    print("dE = ", dE)
#    v0 = v1
#    #iterShoot = iterShoot + 1

#print("Outside E = ", E)
#v1 = f_root_sch(E, l, Ra, Veff)
#print("v1 = ", v1)

dE = abs(E)/fraction
E = E + dE
print("E = ", E)
v1 = f_root_sch(E, l, Ra, Veff)
print("v1 = ", v1)