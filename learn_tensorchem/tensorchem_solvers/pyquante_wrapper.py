import numpy as np
import time
from math import pi
import tucker3d as tuck

#from pyquante2 import molecule, SCF, CGBF
from pyquante2 import molecule, rhf, basisset
import numpy as np
from numpy.linalg import svd, cholesky, solve


def run_pyquante(mol, x, eps, basis = "sto-3g"):

    atoms = [None]*mol.num_atoms
    for i in range(mol.num_atoms):
        atoms[i] = (mol.charge[i], mol.rad[i,0], mol.rad[i,1], mol.rad[i,2])
        
    mol_pq = molecule(atoms, units="Bohr", name=mol.name)
    start = time.time()
    #
    #solver = SCF(mol_pq, method="HF", ConvCriteria=1e-7, MaxIter=40, basis=basis)
    #solver.iterate()
    #
    bfs = basisset(mol_pq, "sto3g")
    solver = rhf(mol_pq, bfs)
    ens = solver.converge()
    end = time.time()
    print("SCF iteration - done")
    print("Time taken for convergence: ", end - start)
    print("ens = ", ens)

    #norb = int(solver.solver.nclosed)
    #norb = len(solver.orbe)
    norb = mol.orbitals

    cf = solver.orbs[:,0:norb]

    psi = [None]*norb
    for i in range(norb):
        psi[i] = gto2tuck(bfs, cf[:,i], x, eps)
        print(psi[i].r)
        #psi[i] = local(psi[i],1e-14)
        #psi[i] = tuck.round(psi[i],eps)
    E = solver.orbe[:norb]
    return psi, E



def local(a,eps):
    
    b = tuck.zeros(a.n, dtype = np.float64)
    b.n = a.n

    u1,v1,r1 = svd_trunc(a.u[0],eps)
    u2,v2,r2 = svd_trunc(a.u[1],eps)
    u3,v3,r3 = svd_trunc(a.u[2],eps)
    
    b.u[0] = np.real(u1)
    b.u[1] = np.real(u2)
    b.u[2] = np.real(u3)
    b.r = (r1,r2,r3)
    

    g = np.dot(a.core, np.transpose(v3))
    g = np.transpose(g, [2,0,1])
    g = np.dot(g, np.transpose(v2))
    g = np.transpose(g, [0,2,1])
    g = np.dot(g, np.transpose(v1))
    b.core = np.real(np.transpose(g, [2,1,0]))
  


    return b



def svd_trunc(A, eps = 1e-14):
    
    u, s, v = np.linalg.svd(A,full_matrices = False)
    
    N1, N2 = A.shape
    
    eps_svd = eps*s[0]/np.sqrt(3)
    r = min(N1, N2)
    for i in range(min(N1, N2)):
        if s[i] <= eps_svd:
            r = i
            break
    #print s/s[0]
    u = u[:,:r].copy()
    v = v[:r,:].copy()
    s = s[:r].copy()
    
    return u, np.dot(np.diag(s),(v)), r


def gto2tuck(w, cf, x, eps):
    
    n = len(x)
    r = len(cf)
    
    k = 0
    for alpha in range(r):
        k += len(w[alpha].pgbfs)
    print("k = ", k)
    
    U1 = np.zeros((n,k*r))
    U2 = np.zeros((n,k*r))
    U3 = np.zeros((n,k*r))
    c = np.zeros(k*r)
    
    print("cf = ", cf)

    nG = 0
    for alpha in range(r):
        for beta in range(len(w[alpha].pgbfs)):
            prim = w[alpha].pgbfs[beta]
            i, j, k = prim.powers
            x0, y0, z0 = prim.origin
            U1[:,beta + nG] = pow(x-x0,i)*np.exp(-prim.exponent*(x-x0)**2)
            U2[:,beta + nG] = pow(x-y0,j)*np.exp(-prim.exponent*(x-y0)**2)
            U3[:,beta + nG] = pow(x-z0,k)*np.exp(-prim.exponent*(x-z0)**2)
            #c[beta + nG] = w[alpha].norm * prim.norm * prim.coef*cf[alpha]
            c[beta + nG] = prim.norm * w[alpha].coefs[beta] * cf[alpha]
        nG += len(w[alpha].pgbfs)

    print("sum U1 = ", U1.sum())
    print("sum U2 = ", U2.sum())
    print("sum U3 = ", U3.sum())

    print("Doing cross2d_full")
    u1, v1 = tuck.cross.cross2d_full(U1, eps)
    u2, v2 = tuck.cross.cross2d_full(U2, eps)
    u3, v3 = tuck.cross.cross2d_full(U3, eps)
    print("Pass here 126 in pyquante wrappers")
            
    v1 = np.real(H(v1))
    v2 = np.real(H(v2))
    v3 = np.real(H(v3))
            
    r1 = v1.shape[0]
    r2 = v2.shape[0]
    r3 = v3.shape[0]

    maxrank = max([r1, r2, r3])
    core = np.zeros((maxrank, maxrank, maxrank))
    for a1 in range(r1):
        for a2 in range(r2):
            for a3 in range(r3):
                for a in range(len(c)):
                    core[a1, a2, a3] += c[a] * v1[a1, a] * v2[a2, a] * v3[a3, a]

    u1_new = np.zeros((n, maxrank))
    u2_new = np.zeros((n, maxrank))
    u3_new = np.zeros((n, maxrank))

    u1_new[:, :r1] = u1
    u2_new[:, :r2] = u2
    u3_new[:, :r3] = u3

    b = tuck.tensor()
    b.core = core.copy()
    b.u[0] = np.real(u1_new)
    b.u[1] = np.real(u2_new)
    b.u[2] = np.real(u3_new)

    b.r = core.shape
    b.n = (n, n, n)

    return b



def H(A):
    return np.transpose(np.conjugate(A))
