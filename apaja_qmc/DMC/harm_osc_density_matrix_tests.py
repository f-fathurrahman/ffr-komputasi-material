#
# Compare the three forms of harmonic oscilaltor density matrix
#
import numpy as np
from numpy import sinh, cosh, tanh

m = 1.0
hbar= 1.0
omega = 2.344
arg = 1.2345

def rho_1(x,xp):
    return -(m*omega)/(2*hbar*sinh(arg)) * ((xp**2 + x**2)*cosh(arg) -2*xp*x)


def rho_2(x,xp):
    return  -(m*omega)/(2*hbar)*tanh(arg/2)*(xp**2 + x**2) - (m*omega)/(2*hbar*sinh(arg))*(xp-x)**2

def rho_3(x,xp):
    return -(1/( (2*hbar)/(m*omega) *tanh(arg))) * (xp-x/cosh(arg))**2 - x**2 * (m*omega)/(2*hbar)*tanh(arg)

xs = np.random.random(10)
xps = np.random.random(10)

for x in xs:
    for xp in xps:
        if np.isclose(rho_1(x,xp),rho_2(x,xp)) and  np.isclose(rho_1(x,xp),rho_3(x,xp)):
            print("rho_1 = rho_2 = rho_3")
        else:
            print("not equal:")
            print(rho_1(x,xp),rho_2(x,xp), rho_3(x,xp))

