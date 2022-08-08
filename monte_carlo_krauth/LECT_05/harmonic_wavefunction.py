from math import sqrt, pi
import numpy as np

Nstates = 4
Nx = 51
x = np.linspace(-5.0, 5.0, Nx)
psi = np.zeros((Nx, Nstates))

psi[:,0] = np.exp( -0.5*x**2)/sqrt(pi)
psi[:,1] = sqrt(2.0) * x * psi[:,0]

# other excited states (through recursion):
for n in range(2,Nstates):
    psi[:,n] = sqrt(2.0/n) * x * psi[:,n-1] - sqrt((n - 1)/n) * psi[:,n-2]
