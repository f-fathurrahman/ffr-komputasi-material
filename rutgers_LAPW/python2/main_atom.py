import numpy as np
import matplotlib.pyplot as plt
from math import exp, pi
import scipy

from excor import ExchangeCorrelation
from Numerov import Numerov
from solve_scheq import CRHS

# Some parameters
RmaxAtom = 10.0  # The end of the radial mesh (maximum r)
Nratom = 3001    # Number of points in radial mesh

Z = 29
# Core specification
core = [3,2,0]  # 1s,2s,3s, 1p,2p, no-d

XC = ExchangeCorrelation(3)

# Radial mesh, equally spaced
R0 = np.linspace(1e-10, RmaxAtom, Nratom) 
# Inverse radial mesh (inverted in sequence)
Ra = R0[::-1]

Veff = -np.ones(len(Ra), dtype=float)/Ra

#plt.clf()
#plt.plot(R0, np.log(-Veff))
#plt.savefig("IMG_atom_Veff.pdf")
#print(Veff)

# We add one more state to core to get atomic states
catm = [c + 1 for c in core]
print("catm = ", catm)



def integ_scheq_numerov(Ex, l, R, Veff):
    "For searching the core bound states"
    rhs = CRHS(Ex, l, R, Veff)
    h = (R[-1]-R[0])/(len(R)-1.)
    u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))
    extraplt = u[-2]*(2+h**2*rhs[-2])-u[-3]
    return u[-1]

# Setting the arguments manually
core = catm
R = Ra
Veff = Veff
Z = Z

fraction = 4.0
coreRho = np.zeros(len(R), dtype=float)
coreE = 0
coreZ = 0

states = []
print("core = ", core)
for l in range(len(core)):
    print("l = ", l)
    n = 0                         # number of states found
    E = -0.5*Z*Z/(l+1)**2-3.      # here we starts to look for zero
    dE = abs(E)/fraction          # the length of the first step 
    decrease = abs(E)/(abs(E)-dE) # the step will decrease to zero. Check the formula!
    v0 = integ_scheq_numerov(E, l, R, Veff)      # starting value
    while (E < 0) and (n < core[l]):      # we need ncore[l] bound states
        E += dE
        v1 = integ_scheq_numerov(E, l, R, Veff)
        if v1*v0 < 0:
            Energy = scipy.optimize.brentq(integ_scheq_numerov, E-dE, E, args=(l, R, Veff))
            # Density
            rhs = CRHS(Energy, l, R, Veff)
            u = Numerov(rhs, (R[-1]-R[0])/(len(R)-1.), R[0]*exp(-R[0]), R[1]*exp(-R[1]))
            drho = u*u
            norm = abs(scipy.integrate.simps(drho, R ))
            drho *= 1./(norm*4*pi*R**2)
            
            coreRho += drho * (2*(2*l+1.))
            coreE   += Energy*(2*(2*l+1.))
            coreZ   += 2*(2*l+1)
            states.append( (n,l,Energy) )
            n += 1
        dE /= decrease
        v0 = v1

print('   Found core states for (n,l)=[',)
for state in states:
    print('(%2d,%2d)' % state[:2],)
print('] E=[',)
for state in states:
    print('%18.10f,' % state[2],)
print(']')
