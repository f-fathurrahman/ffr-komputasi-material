# ----------------------------------------------------------
# Tight-binding model for H2O molecule
# ----------------------------------------------------------

# import the pythtb module
from pythtb import *
import numpy as np

# geometry: bond length and half bond-angle
b = 1.0
angle = 54.0*np.pi/180

# site energies [O(s), O(p), H(s)]
eos = -1.5
eop = -1.2
eh = -1.0

# hoppings [O(s)-H(s), O(p)-H(s)]
ts = -0.4
tp = -0.3

# define frame for defining vectors: 3D Cartesian
lat = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
]

# define coordinates of orbitals: O(s,px,py,pz) ; H(s) ; H(s)
orb = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [b*np.cos(angle),  b*np.sin(angle), 0.0],
    [b*np.cos(angle), -b*np.sin(angle), 0.0]
]

# define model
H2O = TightBindingModel(0, 3, lat, orb)

H2O.set_onsite([eos, eop, eop, eop, eh, eh])
H2O.set_hop(ts, 0, 4)
H2O.set_hop(ts, 0, 5)
H2O.set_hop(tp*np.cos(angle), 1, 4)
H2O.set_hop(tp*np.cos(angle), 1, 5)
H2O.set_hop(tp*np.sin(angle), 2, 4)
H2O.set_hop(-tp*np.sin(angle), 2, 5)

# print model
H2O.display()

# solve model
(evals, evecs) = H2O.solve_all(eig_vectors=True)

# the model is real, so OK to discard imaginary parts of eigenvectors
evecs = evecs.real

# optional: choose overall sign of evec according to some specified rule
# (here, we make the average oxygen p component positive)
for i in range(len(evals)):
    if sum(evecs[i,1:4]) < 0:
        evecs[i,:] = -evecs[i,:]

# print results, setting numpy to format floats as xx.xxx
np.set_printoptions(formatter={'float': '{: 6.3f}'.format})

# print eigenvalues and real parts of eigenvectors, one to a line
print(" n eigval eigvec")
for n in range(6):
    print(" %2i %7.3f " % (n, evals[n]), evecs[n,:])
