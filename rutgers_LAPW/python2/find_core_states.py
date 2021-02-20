import numpy as np
from math import exp, pi
import scipy

from solve_scheq import CRHS
from Numerov import Numerov

def find_core_states(core, R, Veff, Z, fraction=4.):
    "Finds all core states"
    def root(Ex, l, R, Veff):
        "For searching the core bound states"
        rhs = CRHS(Ex, l, R, Veff)
        h = (R[-1]-R[0])/(len(R)-1.)
        u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))
        extraplt = u[-2]*(2+h**2*rhs[-2])-u[-3]
        return u[-1]

    coreRho = np.zeros(len(R), dtype=float)
    coreE = 0
    coreZ = 0

    states=[]
    print("in find_core_states: Z = ", Z)
    for l in range(len(core)):
        print("in find_core_states l = ", l)
        n = 0                           # number of states found
        E = -0.5*Z*Z/(l+1)**2-3.      # here we starts to look for zero
        dE = abs(E)/fraction          # the length of the first step 
        decrease = abs(E)/(abs(E)-dE) # the step will decrease to zero. Check the formula!
        v0 = root(E, l, R, Veff)      # starting value
        while E<0 and n<core[l]:      # we need ncore[l] bound states
            E += dE
            v1 = root(E, l, R, Veff)
            if v1*v0<0:
                Energy = scipy.optimize.brentq(root, E-dE, E, args=(l, R, Veff))
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
            dE/=decrease
            v0 = v1

    print('   Found core states for (n,l)=[',)
    for state in states:
        print('(%2d,%2d)' % state[:2],)
    print('] E=[',)
    for state in states:
        print('%18.10f,' % state[2],)
    print(']')
    
    return coreRho[::-1], coreE, coreZ, states
