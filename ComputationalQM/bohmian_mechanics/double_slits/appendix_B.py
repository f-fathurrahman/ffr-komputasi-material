from math import pi
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def psi_y(y,t=0.):
    Y = 1.0
    s0 = 0.2
    k = 0.1
    x = 0.0
    st = s0*( 1.0 + (1j*t)/(2*s0**2) )
    N = (2*pi*st**2)**-0.25
    return N*np.exp( (-((y-Y)**2)/(4*s0*st) ) + 1j*(k*x-((t*k**2)/2)) )

t = np.arange(-0.0, 1.1, 0.01)
y = np.arange(-13.0, 13.0, 0.01)
Y,T = np.meshgrid(y, t)
Z = psi_y(Y,T) + psi_y(-Y,T)

fig = plt.figure(figsize=(17,10))
plt.suptitle("Imaginary surface over time of the wave function for a double slit system", fontsize=20)
ax = fig.gca(projection="3d")
ax.set_zlim([-1.0,1.0])

surf = ax.plot_surface(Y, T, np.real(Z), cmap=matplotlib.cm.jet, linewidth=0, antialiased=False)
ax.set_xlabel("y",fontsize=20)
ax.set_ylabel("t",fontsize=20)
ax.set_zlabel(r"$\psi$",fontsize=30)
ax.view_init(elev=36, azim=63)

cbar = plt.colorbar(surf)
cbar.set_label("Surface gradient", fontsize=20)

plt.savefig("TEMP_appendix_B.pdf")