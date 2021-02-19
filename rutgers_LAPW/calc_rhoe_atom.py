import numpy as np
import scipy
from math import sqrt, exp, pi

from my_utilities import extrapolate
from Numerov import Numerov
from solve_scheq import CRHS

def calc_rhoe_atom(states, R, Veff, Z):
    "Computes electron charge density, given the bound states and Z"
    rho = np.zeros(len(R), dtype=float)
    N = 0
    Ebs = 0
    for state in states:
        l = state[1]
        E = state[2]
        rhs = CRHS(E, l, R, Veff)
        h = (R[-1]-R[0])/(len(R)-1.)
        u = Numerov(rhs, h, R[0]*exp(-R[0]), R[1]*exp(-R[1]))
        u /= sqrt(abs(sum(u)*h))   # To avoid overflow
        
        u2 = u**2
        norm = abs(scipy.integrate.simps(u2, R))
        u2 *= 1./norm

        dN = 2*(2*l+1.)
        if N+dN<Z:
            ferm = 1.
        else:
            ferm = (Z-N)/dN
        drho = u2[:]*dN*ferm/(4*pi*R**2)
        
        rho += drho

        N += dN
        Ebs += E * ferm*dN
        if N>=Z: break

    rho[-1] = extrapolate(0.0, R[-2], R[-3], rho[-2], rho[-3])
    return (rho, Ebs)