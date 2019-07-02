from math import pi
import numpy as np
import matplotlib.pyplot as plt

def psi_y(y, t = 1.0):
    Y = 1.0
    s0 = 0.2
    k = 0.1
    x = 0.0
    st = s0*(1 + (1j*t)/(2*s0**2))
    N = (2*pi*st**2)**0.25
    #
    return N*np.exp((-((y-Y)**2)/(4*s0*st)) + 1j*(k*x-((t*k**2)/2)))


y = np.arange(-13.0, 13.0, 0.01)
Z = psi_y(y) + psi_y(-y)

plt.figure(figsize=(12,10))

plt.subplot(211)
plt.xlabel("y",fontsize=20)
plt.ylabel(r"$\psi$",fontsize=30)
plt.plot(y, np.real(Z), "b", label="Real part")
plt.legend()

plt.subplot(212)
plt.xlabel("y",fontsize=20)
plt.ylabel(r"$\psi$",fontsize=30)
plt.plot(y, np.imag(Z), "r", label="Imaginary part")
plt.legend()

plt.savefig("TEMP_appendix_A.pdf")