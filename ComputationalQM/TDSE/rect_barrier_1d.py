# File: ./examples1D/Rectangular_Barrier_1D.py
# Run as    python3 bpm.py Rectangular_Barrier_1D 1D
# A Gaussian wave packet impinges on a rectangular barrier. Part of the wave
# is reflected, part of the wave is transmitted.

import numpy as np

# Grid points
Nx = 1600 
Ny = Nx

# Evolution step
dt = 0.0001

# Propagation end
tmax = 6

# x- and y-window size
xmax = 50
ymax = xmax

# number of .png images
images = 100

# 0 = periodic boundary
absorb_coeff = 20

# If 1, it plots on the screen but does not save the images
# If 2, it saves the images but does not plot on the screen
# If 3, it saves the images and plots on the screen
output_choice = 3

# Fixes a maximum scale of |psi|**2 for the plots. If 0, it does not fix it.
fixmaximum= 1.05

# Initial wavefunction
# A Gaussian wave packet moving rightwards
def psi_0(x,y):
    vx = 10     # value of the initial velocity
    f = 0.j+np.exp(-((x+15)/4)**2)  # Gaussian profile
    f = f*np.exp(1.j*vx*x)   # Multiply by an x-dependent phase to introduce velocity
    return f

# A barrier modeled by V=0 for |x|>5 and V=40 for |x|<5
def V(x,y,t,psi):
    V = np.piecewise(x, [abs(x-5)<2.5, abs(x-5)>=2.5],[40,0])
    return V;

