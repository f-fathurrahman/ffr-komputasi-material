#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
from solvers import solver
import numpy as np
import qtools as qt
import tucker3d as tuck


# In[17]:


name = 'H2'
molecule = qt.molecule(name)

gradients_005 = []
gradients_01 = []
gradients_02 = []
gradients_03 = []
gradients_05 = []
sols = []
for n in [2**11, 2**12, 2**13, 2**14]:
    sol = solver(molecule, method='hf', mixing=2, eps=1e-11,
                 maxiter=60, meshsize=n, boxsize=12.0)
    
    sol.solve()
    gradients_005.append(gradient_diatomic(sol, smoothing=0.005))
    gradients_01.append(gradient_diatomic(sol, smoothing=0.01))
    gradients_02.append(gradient_diatomic(sol, smoothing=0.02))
    gradients_03.append(gradient_diatomic(sol, smoothing=0.03))
    gradients_05.append(gradient_diatomic(sol, smoothing=0.05))
    
    print gradients


# In[18]:


print gradients_005
print gradients_01
print gradients_02
print gradients_03
print gradients_05


# In[19]:


def aitken(a2, a1, a0):
    return a2 - (a2-a1)**2/(a2-2*a1+a0)

print aitken(-0.005391281864745, -0.005391911296848, -0.0054030852009102)


# In[16]:


print gradients_005
print gradients_01
print gradients_02
print gradients_03
print gradients_05


# In[9]:


gradient_diatomic(sol, smoothing=0.01)


# In[12]:


gradient_diatomic(sol, smoothing=0.005)


# In[33]:


from PyQuante.Molecule import Molecule
from PyQuante.dft import dft
h2 = Molecule('h2',[(1,(0,0,0)),(1,(1.4,0,0))])


# In[35]:


en,orbe,orbs = dft(h2,functional='LDA')


# In[ ]:





# In[8]:


gradient_diatomic(sol, smoothing=0.2)


# In[14]:


2*a - b


# In[4]:


2*a - b


# In[16]:


x = sol.params.grid


# In[17]:


tuck.dot(sol.rho, tuck.ones(sol.rho.n))*(x[1]-x[0])**3


# In[4]:


from scipy.special import erf

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
        def fun_coulomb_squared((i, j, k)):
            r = np.sqrt((vec[0] - x[i])**2 + (vec[1] - x[j])**2 + (vec[2] - x[k])**2)
            return charge*(x[i] - vec[0]) / r * (-smoothed_coulomb_derivative(r/smoothing))/smoothing**2
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
    
def smoothed_coulomb_derivative(r):
        if r>=6:
            return -1.0/r**2
        elif r<6 and r>=0.1:
            return 2*np.exp(-r**2)/np.sqrt(np.pi)/r - erf(r)/r**2 - 1.0/3/np.sqrt(np.pi) * (2*r*np.exp(-r**2) + 128*r*np.exp(-4*r**2))
        else:
            return -4.0/3*r + 4.0/5*r**3 - 2.0/7*r**5 + 2.0/27*r**7 - 1.0/66*r**9 - 1.0/3/np.sqrt(np.pi) * (2*r*np.exp(-r**2) + 128*r*np.exp(-4*r**2))


# In[17]:


from scipy.special import erf

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
        def fun_coulomb_squared((i, j, k)):
            r = np.sqrt((vec[0] - x[i])**2 + (vec[1] - x[j])**2 + (vec[2] - x[k])**2)
            return charge*(x[i] - vec[0]) / r * (-smoothed_coulomb_derivative(r/smoothing))/smoothing**2
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


        return tuck.dot(d_coulomb, 2*tuck.ones(d_coulomb.n))*(x[1]-x[0])**3 #+ dE_nuc
    
def smoothed_coulomb_derivative(r):
        if r>=6:
            return -1.0/r**2
        elif r<6 and r>=0.1:
            return 2*np.exp(-r**2)/np.sqrt(np.pi)/r - erf(r)/r**2 - 1.0/3/np.sqrt(np.pi) * (2*r*np.exp(-r**2) + 128*r*np.exp(-4*r**2))
        else:
            return -4.0/3*r + 4.0/5*r**3 - 2.0/7*r**5 + 2.0/27*r**7 - 1.0/66*r**9 - 1.0/3/np.sqrt(np.pi) * (2*r*np.exp(-r**2) + 128*r*np.exp(-4*r**2))


# In[10]:


def gradient_diatomic(self):

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
        def fun_coulomb_squared((i, j, k)):
            return charge*(x[i] - vec[0]) / np.sqrt((vec[0] - x[i])**2 + (vec[1] - x[j])**2 + (vec[2] - x[k])**2)**(3.0)
        d_coulomb = tuck.cross.cross3d(fun_coulomb_squared, N, self.params.eps)
        d_coulomb = tuck.round(tuck.real(d_coulomb), self.params.eps)

        self.get_rho()
        #d_coulomb = tuck.cross.multifun([self.rho, d_coulomb], self.params.eps, lambda x: x[0]*x[1])
        
        
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

