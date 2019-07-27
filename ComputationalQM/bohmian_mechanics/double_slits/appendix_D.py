from math import pi, sqrt
import numpy as np
import matplotlib.pyplot as plt

def wave(x, y, t=15.0, a=36.0, c=1.0, theta=pi/4.0):
    k0 = sqrt(2.0)*pi/2.0
    sint = np.sin(theta)
    cost = np.cos(theta)
    return ( np.exp( 1.0j*(k0*(x*cost - y*sint) - c*t) ) *
             np.exp( -( (x*cost - y*sint)-c*t)**2/a ) *
             np.exp( -(x*sint + y*cost)**2/a ) +
             np.exp( 1.0j*(k0*(x*cost+y*sint) - c*t) ) *
             np.exp( -( (x*cost + y*sint) - c*t)**2/a ) *
             np.exp( -(x*sint - y*cost)**2/a) )

delta = 0.1

x = np.arange(-12.0, 12.0, delta)
y = np.arange(-12.0, 12.0, delta)

X, Y = np.meshgrid(x, y)

plt.figure(figsize=(14,4))
plt.suptitle('Contour plot of two interfering wave packets', fontsize=15)

ax = plt.subplot(131)

plt.text(0.9, 0.9,"t=-15s", ha="center", va="center", transform=ax.transAxes)
plt.xlabel("x")
plt.ylabel("y")

ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])

Z = wave(X,Y,t=-15.)
plt.contour( np.real(Z*np.conj(Z)) )

ax = plt.subplot(132)
plt.text(0.9, 0.9,"t=0s", ha="center", va="center", transform=ax.transAxes)
plt.xlabel("x")
plt.ylabel("y")
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
Z = wave(X, Y, t=0)
plt.contour( np.real(Z*np.conj(Z)) )


ax = plt.subplot(133)
plt.text(0.2, 0.9,"t=15s", ha="center", va="center", transform=ax.transAxes)
plt.xlabel("x")
plt.ylabel("y")
ax.xaxis.set_ticks([])
ax.yaxis.set_ticks([])
Z = wave(X,Y,t=15.)
c = plt.contour( np.real(Z*np.conj(Z)) )

cbar = plt.colorbar(c)
cbar.set_ticks([])

plt.savefig("TEMP_appendix_D.pdf")
