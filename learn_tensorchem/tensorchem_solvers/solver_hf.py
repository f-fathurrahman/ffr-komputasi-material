import numpy as np
import tucker3d as tuck
from qtools import pots
from qtools.basic import bilinear as qt_bilinear
from qtools.basic import UT_prod as qt_UT_prod
from qtools.basic import prod as qt_prod
from qtools.basic import normalize as qt_normalize
from qtools.poisson import poisson as qt_poisson


# m is number of mixed vectors

def solver_hf(solver, molecule, psi, E, eps, grid, T, ind, max_iter, pot_coulomb, poisson_solver="Fourier"):
    
    a = grid[0]
    b = grid[1]
    N = grid[2]
    
    h = (b-a)/(N-1)
    
    x = np.zeros(N)
    for i in range(N):
        x[i] = a + i*h

    solver.iterative.orb_energies = []
    solver.iterative.convergence = []
    
    k=0
    for i in range(max_iter):
        psi, E_new = hf1iter(molecule, psi, E, grid, T, pot_coulomb, eps, poisson_solver)
        #print E
        err = abs((E-E_new)/E)
        solver.iterative.orb_energies.append(E_new)
        solver.iterative.convergence.append(err)
        print(f"Iteration = {i} accuracy = {max(err):.8g} eps={eps}")
        if max(err) < eps:
            print("Process converged with", i, "iterations")
            break
            #k += 1
            #if k == 4:
            #    print("Process converged with", i, "iterations")
            #    break
        
        E = E_new
        #print E

    if i == max_iter - 1:
        print("")
        print("Process did not converge with eps precision")
        print("")

    solver.orb_energies = E_new
    E_full = hf_full_energy(molecule, psi, E, grid, T, eps)
    solver.energy = E_full
    solver.psi = psi
    print("Full energy = ", E_full)
    #return psi, E, E_full



def hf1iter(molecule, psi, E, grid, T, pot_coulomb, eps, poisson_solver, num_iter=1):
    
    
    E_correct = molecule.energy
    
    eps_exchange = eps
    
    a = grid[0]
    b = grid[1]
    N = grid[2]
    
    h = (b-a)/(N-1)
    
    x = np.zeros(N)
    for i in range(N):
        x[i] = a + i*h
    
    
    sf = h ** (3./2) # scaling factor
    
    Norb = molecule.orbitals
    #Norb = len(psi)
    print("Norb = ", Norb)

    V = [0]*Norb
    psi_new = [0]*Norb


    prod = lambda A,B,C: tuck.cross.multifun([A,B], eps, lambda a_b: a_b[0]*a_b[1], y0=C)
    
    for k in range(num_iter):

        
        density = tuck.zeros((N,N,N))  # density calculation
        for i in range(Norb):
            density = density + prod(tuck.conj(psi[i]), psi[i], psi[i]) # !can be faster!
            density = tuck.round(density, eps)
        pot_hartree = tuck.cross.conv(T, density, eps,pr=None)
        pot_hartree = tuck.round(tuck.real(pot_hartree),eps)

    
        E_electr = 0
        for i in range(Norb):
            V[i] =  prod(tuck.round(pot_coulomb + 2*pot_hartree, eps), psi[i], psi[i])
            exchange = pots.exchange(psi, psi, i, eps_exchange, T, molecule)
            V[i] = tuck.round(V[i] - exchange, eps)
            #
            # Full energy computation
            #V_energ = tuck.round(prod(2*pot_hartree, psi[i],psi[i]) - exchange, eps)
            #E_electr += 2*E[i] - sf**2*tuck.dot(psi[i], V_energ)
            psi_new[i] = -qt_poisson(2*V[i], -2*E[i], N, h, eps, solver=poisson_solver)
            #if abs( (E_full_1-E_full_0)/E_full_1 ) <= 3 * eps:
            #check += 1
                #if check == 3:
        #break

        
        ################# Fock matrix ###################
        res1 = qt_bilinear(psi_new, psi_new)*sf**2
        print("res1 = ", res1)
        L = np.linalg.cholesky(res1) # orthogonalization
        psi_Q = qt_UT_prod(psi_new, H(np.linalg.inv(L)), eps)
        
        # Fock matrix
        
        V_new = [0]*Norb
        for i in range(Norb):
            V_new[i] = prod(tuck.round(pot_coulomb + 2*pot_hartree,eps), psi_Q[i],psi_Q[i])
            exchange_new = pots.exchange(psi_Q, psi, i, eps_exchange, T, molecule)
            V_new[i] = tuck.round(V_new[i] - exchange_new, eps)

        print("E = ", E)
        print("L = ", L)

        Energ1 = np.dot(np.diag(E), H(np.linalg.inv(L)))

        Energ2 = qt_bilinear(psi_new, psi_new)*sf**2
        Energ2 = np.dot(np.linalg.inv(L), Energ2)
        Energ = np.dot(Energ2, Energ1)
        F = qt_bilinear(psi_Q, V_new)*sf**2 - np.dot(qt_bilinear(psi_Q, V)*sf**2, H(np.linalg.inv(L))) + Energ
        #F = np.real((F + F.T)/2)
        
        E_new = np.zeros(Norb, dtype = np.float64)
        E_new, S = np.linalg.eigh(F)
        
        
        
        psi = qt_prod(psi_Q, S.T, density, eps)
        for i in range(Norb):
            psi[i] = qt_normalize(psi[i], h)
            psi[i] = tuck.round(psi[i], eps)
        
        E = E_new.copy()
        
        for i in range(Norb):
            if E[i]>0:
                E[i] = -0.0
    
    
    return psi, E



def hf_full_energy(molecule, psi, E, grid, T, eps):
    
    eps_exchange = eps
    
    a = grid[0]
    b = grid[1]
    N = grid[2]
    
    h = (b-a)/(N-1)
    
    sf = h ** (3./2)
    
    Norb = molecule.orbitals
    num_atoms = molecule.num_atoms
    
    prod = lambda A,B,C: tuck.cross.multifun([A,B], eps, lambda a_b1: a_b1[0]*a_b1[1], y0=C)
    
    density = tuck.zeros((N,N,N))  # density calculation
    for i in range(Norb):
        density = density + prod(psi[i], psi[i], psi[i]) # !can be faster!
        density = tuck.round(density, eps)
    #
    pot_hartree = tuck.cross.conv(T, density, eps)
    pot_hartree = tuck.real(pot_hartree)
    pot_hartree = tuck.round(pot_hartree, eps)
    #
    V = [0] * Norb
    for i in range(Norb):
        V[i] =  (prod(tuck.round(2*pot_hartree, eps), psi[i], psi[i]) -
                 pots.exchange(psi, psi, i, eps_exchange, T, molecule))
        V[i] = tuck.round(V[i], eps)
    #
    E_electr = 0
    for i in range(Norb):
        E_electr = E_electr + 2*E[i] - sf**2*tuck.dot(psi[i], V[i]) 
    #
    E_nuc = 0
    for i in range(num_atoms):
        vec_i = molecule.rad[i]
        charge_i = molecule.charge[i]
        for j in range(i):
            vec_j = molecule.rad[j]
            charge_j = molecule.charge[j]
            E_nuc = (E_nuc + charge_i*charge_j/
                     np.sqrt((vec_i[0] - vec_j[0])**2 + (vec_i[1] - vec_j[1])**2 + (vec_i[2] - vec_j[2])**2))
    return E_nuc + E_electr


def H(A):
    return np.transpose(np.conjugate(A))


