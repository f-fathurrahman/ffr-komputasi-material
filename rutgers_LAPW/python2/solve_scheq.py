import scipy
import weave
import numpy as np
from math import sqrt
from scipy import integrate

from my_utilities import extrapolate
from Numerov import Numerov, NumerovGen
from CRHS import CRHS

def solve_scheq(Z, Enu, R0, RMuffinTin, Veff):
    """Solves the SCH Eq for Psi(Enu) and its energy derivative
       Returns logarithmic derivative, Psi(l,E) and its energy derivative
       
       Please  see Eq.30 on page 20 for definition, and Eq.49 on page 26
        for S*Psi'(S)/Psi(S) and S*dPsi'(S)/dPsi(S)
    """
    def startSol(Z, l, r):
        "good choice for starting Numerov algorithm"
        return r**(l+1)*(1-Z*r/(l+1))

    logDer=[]
    Psi_l=[]
    Psip_l=[]
    for l in range(len(Enu)):

        # Computes Psi=u/r
        crhs = CRHS(Enu[l], l, R0, Veff)
        crhs[0]=0
        ur = Numerov(crhs, (R0[-1]-R0[0])/(len(R0)-1.), 0.0, startSol(Z,l,R0[1]))
        
        ur *= 1/sqrt( scipy.integrate.simps(ur*ur, R0) )  # normalization
        Psi_l.append( ur/R0 ) # storing Psi
        Psi_l[-1][0] = extrapolate(R0[0], R0[1], R0[2], ur[1]/R0[1], ur[2]/R0[2])
        
        # For energy derivative of Psi' = urp/r
        inhom = -2*ur
        urp = NumerovGen(crhs, inhom, (R0[-1]-R0[0])/(len(R0)-1.), 0.0, startSol(Z,l,R0[1]))

        # Energy derivative should be orthogonal
        alpha = scipy.integrate.simps(ur*urp, R0)
        urp -= alpha*ur
        Psip_l.append( urp/R0 ) # storing Psip'
        Psip_l[-1][0] = extrapolate(R0[0], R0[1], R0[2], urp[1]/R0[1], urp[2]/R0[2])
        
        # <\cdot{\psi}|\cdot{\psi}>
        PsipPsip = integrate.simps(urp*urp, R0)
        
        # Computes the logarithmic derivative
        v1 = crhs[-1]*ur[-1]
        v0 = crhs[-2]*ur[-2]
        w1 = crhs[-1]*urp[-1]+inhom[-1]
        w0 = crhs[-2]*urp[-2]+inhom[-2]
        dh = R0[2]-R0[1]
        dudr  = (ur[-1]-ur[-2])/dh + 0.125*dh*(3*v1+v0)
        dupdr = (urp[-1]-urp[-2])/dh + 0.125*dh*(3*w1+w0)
        
        dlogPsi = RMuffinTin*dudr/ur[-1] - 1
        dlogPsip = RMuffinTin*dupdr/urp[-1] - 1
        Psi = ur[-1]/RMuffinTin
        Psip = urp[-1]/RMuffinTin
        
        logDer.append( (Psi, Psip, dlogPsi, dlogPsip, PsipPsip) )
        
    return (logDer, Psi_l, Psip_l)

