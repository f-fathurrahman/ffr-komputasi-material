import numpy as np
from math import sqrt
from numpy.linalg import eigh

def do_optim_BFGS(atoms, fmax_conv=0.01):

    # Parameter
    maxstep = 0.04
    
    # Initial Hessian
    H = np.eye(3 * len(atoms)) * 70.0
    #
    r0 = None
    f0 = None

    for iterBFGS in range(0,100):

        print()
        print("Start iterBFGS = ", iterBFGS)
        print()

        Etot = atoms.get_potential_energy()
        forces = atoms.get_forces()
        r = atoms.get_positions()
        f = forces.reshape(-1) # linearize

        print("Etot = ", Etot)
        print("f (linear index) = ")
        print(f)

        fmax = sqrt((forces ** 2).sum(axis=1).max())
        print("fmax = ", fmax)
        if fmax < fmax_conv:
            print()
            print("*** BFGS convergence is achieved ***")
            print()
            break

        if iterBFGS > 0:
            H = bfgs_update(H, r.flat, f, r0, f0)
        print("H = ")
        print(H)

        omega, V = eigh(H)
        print("omega = ")
        print(omega)
        
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        print("dr = ")
        print(dr)
    
        steplengths = (dr**2).sum(1)**0.5
        dr = bfgs_determine_step(maxstep, dr, steplengths)
        print("dr after bfgs_determine_step= ")
        print(dr)
    
        print("Update positions")
        atoms.set_positions(r + dr)
    
        r0 = r.flat.copy()
        f0 = f.copy()

    return

def bfgs_update(H, r, f, r0, f0):
    dr = r - r0
    if np.abs(dr).max() < 1e-7:
        # Same configuration again (maybe a restart):
        return
    df = f - f0
    a = np.dot(dr, df)
    dg = np.dot(H, dr)
    b = np.dot(dr, dg)
    H = H - (np.outer(df, df) / a + np.outer(dg, dg) / b)
    return H


# Determine step to take according to maxstep
# Normalize all steps as the largest step. This way
# we still move along the eigendirection.
def bfgs_determine_step(maxstep, dr, steplengths):
    maxsteplength = np.max(steplengths)
    if maxsteplength >= maxstep:
        dr = dr * maxstep / maxsteplength
    return dr
