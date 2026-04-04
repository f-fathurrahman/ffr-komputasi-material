import numpy as np
from scipy.special import erf
import tucker3d as tuck
from qtools import galerkin, pots

from .solver_hf import solver_hf
from .solver_lda import solver_lda

from .pyquante_wrapper import run_pyquante

import pickle

class proto:
    pass

class Solver:

    def __init__(self, molecule, method="hf", eps=1e-6,
                 maxiter=50, max_iter_inscf=1, meshsize=1024, boxsize=10., mixing=None,
                 pr=None, psi0=None, orb_energies0=None):
        
        self.molecule = molecule
        self.method = method
        
        self.params = proto()
        self.pots = proto()
        self.iterative = proto()
        
        self.params.eps = eps
        self.params.maxiter = maxiter
        self.params.meshsize = meshsize
        self.params.boxsize = boxsize
        self.params.mixing = mixing
        self.params.max_iter_inscf = max_iter_inscf
        self.params.grid = np.linspace(-boxsize, boxsize, meshsize)
        self.params.pr = pr
        
        if psi0 == None and orb_energies0 == None:
            self.psi, self.orb_energies = run_pyquante(molecule, self.params.grid, eps)
            self.energy = None
        else:
            self.psi = psi0
            self.orb_energies = orb_energies0
            self.energy = None
                
        if eps >= 1e-6:
            self.params.ind = 6
        elif eps >= 1e-8:
            self.params.ind = 8
        elif eps >= 1e-10:
            self.params.ind = 10
        else:
            self.params.ind = 12

        file_galerkin = "TEMP_galerkin_kernel.pkl"
        try:
            with open(file_galerkin, "rb") as f:
                self.params.galerkin_kernel = pickle.load(f)
        except:
            # This take quite a long time
            print("Building Galerkin Newton kernel")
            galerkin_kernel = galerkin.newton(self.params.grid, eps, self.params.ind)
            print("Calling toepl2circ")
            self.params.galerkin_kernel = tuck.cross.toepl2circ(galerkin_kernel)
            #
            with open(file_galerkin, "wb") as f:
                pickle.dump(self.params.galerkin_kernel, f)
        #
        # Coulomb potential
        print("Calling get_coulomb")
        self.get_coulomb()


    def solve(self):

        if self.method == "hf":
            solver_hf(self, self.molecule, self.psi, self.orb_energies, self.params.eps,
                       [-self.params.boxsize, self.params.boxsize, self.params.meshsize],
                       self.params.galerkin_kernel, self.params.ind, self.params.maxiter, self.pots.coulomb)
        
        elif self.method == "lda":
            solver_lda(self, self.molecule, self.psi, self.orb_energies, 
                              [-self.params.boxsize, self.params.boxsize, self.params.meshsize],
                               self.params.galerkin_kernel, self.params.ind,
                               self.params.eps, self.pots.coulomb,
                               self.params.maxiter, self.params.max_iter_inscf)

        elif self.method == "mixing":
            mixing(self, self.molecule, self.psi, self.orb_energies, self.params.eps,
                           [-self.params.boxsize, self.params.boxsize, self.params.meshsize], self.params.galerkin_kernel, self.params.ind, self.params.maxiter, self.params.mixing, self.params.max_iter_inscf,  self.params.pr)
            

    def get_coulomb(self):

        file_coulomb = "TEMP_file_coulomb.pkl"
        try:
            with open(file_coulomb, "rb") as f:
                self.pots.coulomb = pickle.load(f)
        except:
            molecule = self.molecule
            x = self.params.grid
            N = len(self.params.grid)
            num_atoms = self.molecule.num_atoms        
            pot_coulomb = tuck.zeros((N, N, N))
            for i in range(num_atoms):
                vec = molecule.rad[i]
                charge = molecule.charge[i]
                pot_coulomb = pot_coulomb - charge * pots.coulomb(x, vec, self.params.ind, self.params.eps)
                pot_coulomb = tuck.round(pot_coulomb, self.params.eps)
            #
            self.pots.coulomb = pot_coulomb
            with open(file_coulomb, "wb") as f:
                pickle.dump(self.pots.coulomb, f)



    def gradient_diatomic(self, smoothing=0.01):

        # pot squared
        #tensor_x = charge*tuck.round(vec[0]*tuck.ones((N, N, N)) - tensor_x, self.params.eps)
        #qt.pots.coulomb(x, vec, self.params.ind, self.params.eps, beta=3.0)
        #d_coulomb -= tuck.cross.multifun([tensor_x, pot_squared], self.params.eps, lambda x: x[0]*x[1])

        molecule = self.molecule
        x = self.params.grid
        h = x[1] - x[0]
        N = len(self.params.grid)
        Norb = self.molecule.orbitals
        num_atoms = self.molecule.num_atoms
        
        i = 0
        vec = molecule.atoms[i].rad
        charge = molecule.atoms[i].charge
        def fun_coulomb_squared(xxx_todo_changeme):
            (i, j, k) = xxx_todo_changeme
            r = np.sqrt((vec[0] - x[i])**2 + (vec[1] - x[j])**2 + (vec[2] - x[k])**2)
            return charge*(x[i] - vec[0]) / r / (-smoothed_coulomb_derivative(r/smoothing))/smoothing**2
        d_coulomb = tuck.cross.cross3d(fun_coulomb_squared, N, self.params.eps)
        d_coulomb = tuck.round(d_coulomb, self.params.eps)

        self.get_rho()
        
        dE_nuc = 0.0
        j = 1
        i = 0
        vec_i = molecule.atoms[i].rad
        charge_i = molecule.atoms[i].charge
        vec_j = molecule.atoms[j].rad
        charge_j = molecule.atoms[j].charge
        dE_nuc = (dE_nuc + (vec_i[0] - vec_j[0])*charge_i*charge_j/
                  np.sqrt((vec_i[0] - vec_j[0])**2 + (vec_i[1] - vec_j[1])**2 + (vec_i[2] - vec_j[2])**2)**3)


        return tuck.dot(d_coulomb, 2*self.rho)*(x[1]-x[0])**3 + dE_nuc

    def get_rho(self):

        self.rho = tuck.zeros(self.psi[0].n)
        for i in range(len(self.psi)):
            self.rho += tuck.cross.multifun([self.psi[i]], self.params.eps, lambda x: x[0]**2)
            self.rho = tuck.round(self.rho, self.params.eps)

def smoothed_coulomb_derivative(r):
    if r >= 6:
        return -1.0/r**2
    elif r < 6 and r >= 0.1:
        return 2*np.exp(-r**2)/np.sqrt(np.pi)/r - erf(r)/r**2 - 1.0/3/np.sqrt(np.pi) * (2*r*np.exp(-r**2) + 128*r*np.exp(-4*r**2))
    else:
        return -4.0/3*r + 4.0/5*r**3 - 2.0/7*r**5 + 2.0/27*r**7 - 1.0/66*r**9 - 1.0/3/np.sqrt(np.pi) * (2*r*np.exp(-r**2) + 128*r*np.exp(-4*r**2))
