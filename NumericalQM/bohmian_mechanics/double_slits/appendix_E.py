import os
import time
from math import sqrt, pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt

def wave(x,y,t=15.):
    k0 = sqrt(2.)*pi/2.0
    a = 3.0
    c = 1.0
    sint = sin(pi/4.)
    cost = cos(pi/4.)
    return ( np.exp(1.j*(k0*(x*cost - y*sint) - c*t))*5.0*sqrt(pi) *
             np.exp(-((x*cost - y*sint)-c*t)**2./4./a**2.0)/a*5.0*sqrt(pi) *
             np.exp(-(x*sint + y*cost)**2./4./a**2.)/a +
             np.exp(1.j*(k0*( x*cost + y*sint)-c*t) )*5.0*sqrt(pi) *
             np.exp(-((x*cost + y*sint)-c*t)**2.0/4.0 /a**2.0)/a*5.*sqrt(pi) *
             np.exp(-(x*sint - y*cost)**2.0/4.0/a**2.)/a )

delta = 0.1
x = np.arange(-12.0, 12.0, delta)
y = np.arange(-12.0, 12.0, delta)
X, Y = np.meshgrid(x, y)

t = -25.0
i = 0
while t < 25.0:
    i = i + 1
    plt.clf()
    Z = wave(X, Y, t)
    plt.contour( np.real(Z*np.conj(Z)) )
    plt.savefig("TEMP_appendix_E_" + str(i) + ".png", dpi=150)
    t = t + 1.0
    print("t = ", t)

