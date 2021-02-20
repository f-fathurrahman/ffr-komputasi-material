import numpy as np
from excor import ExchangeCorrelation
from find_core_states import *
from calc_rhoe_atom import *
from atom_cmpb import *
from solve_poisson import *
from my_utilities import calc_rs

def calc_atom_charge(Z, core, mix=0.3, RmaxAtom=10., Natom=3001, precision=1e-5, Nitt=100):
    """ Computes Atomic electronic density and atomic Energy
    Input:
       Z             --  Nucleolus charge
       core          --  States treated as core in LAPW (example: [3,2,0]  # 1s,2s,3s, 1p,2p, no-d)
       mix           --  Mixing parameter for density
       RmaxAtom      --  The end of the radial mesh (maximum r)
       Natom         --  Number of points in radial mesh
       precision     --  How precise total energy we need      
       Nitt          --  Maximum number of itterations
       """

    XC = ExchangeCorrelation(3)  # VWN

    R0 = np.linspace(1e-10, RmaxAtom, Natom) # Radial mesh
    Ra = R0[::-1]                         # Inverse radial mesh

    Veff = -np.ones(len(Ra), dtype=float)/Ra

    catm = [c + 1 for c in core]          # We add one more state to core to get atomic states
   
    Etot_old = 0
   
    # Finds bound states
    coreRho, coreE, coreZ, states = find_core_states(catm, Ra, Veff, Z)
    
    #exit()

    # Sorts them according to energy
    states.sort(atom_cmpb)

    # Computes charge
    (rho, Ebs) = calc_rhoe_atom(states, Ra, Veff, Z)
    rho = rho[::-1]

    for itt in range(Nitt):

        # Here we have increasing R ->

        # Hartree potential
        UHartree = solve_poisson(Z, R0, rho)
        # Adding exchange-correlation part
        Vxc = [XC.Vx(rsi)+XC.Vc(rsi) for rsi in calc_rs(rho)]
        ExcVxc = [XC.EcVc(rsi) + XC.ExVx(rsi) for rsi in calc_rs(rho)]
        Veff = (UHartree - Z)/R0 + Vxc
        Veff=Veff[::-1]

        # Here we have decreasing R <-
       
        # Finds bound states
        (coreRho, coreE, coreZ, states) = find_core_states(catm, Ra, Veff, Z)

        if itt == 1: exit()  # ffr

        # Sorts them according to energy
        states.sort(atom_cmpb)
        # Computes charge
        (nrho, Ebs) = calc_rhoe_atom(states, Ra, Veff, Z)

        # Total energy
        pot = (ExcVxc*R0**2-0.5*UHartree*R0)*nrho[::-1]*4*pi

        Etot = scipy.integrate.simps(pot, R0) + Ebs
        Ediff = abs(Etot-Etot_old)
       
        print('   %d) Etot=%f Eband=%f Ediff=%f' % (itt, Etot, Ebs, Ediff))
       
        # Mixing
        rho = mix*nrho[::-1] + (1-mix)*rho
        Etot_old = Etot

        if Ediff < precision: break

    return (R0, rho)
