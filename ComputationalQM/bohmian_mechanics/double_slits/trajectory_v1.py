import time
from math import pi
import numpy as np
import sympy
import matplotlib.pyplot as plt

Nparticles = 200
Y = 1.0
sigma_0 = 0.2
kx = 0.1
dt = 1e-2

y = sympy.var("y")
t = sympy.var("t")
x = y*t

st = sigma_0*( 1.0 + (1.0j*t)/(2*sigma_0**2) )
N = (2*pi*st**2)**-0.25

wavefunc_analytic = N*sympy.exp((-((y - Y)**2.0)/(4.0*sigma_0*st)) + 1.0j*(kx*x - ((t*kx**2.0)/2.0)))

print(wavefunc_analytic)

ydot_analytic = wavefunc_analytic.diff(y)/wavefunc_analytic
print(sympy.simplify(ydot_analytic))

I = 1.0j

def ydot(t, y):
    #yy = (Y - y)/(2*sigma_0**2 + 1.0j*t)
    yy = (0.1*I*t*(2.0*I*t + 0.16) - 2.0*(y - 1.0)**1.0)/(2.0*I*t + 0.16)
    return np.imag(yy)

def xdot():
    return kx

def rk4(f, t, y, dt):
    k1 = dt * f(t, y)
    k2 = dt * f(t + 0.5*dt, y + 0.5*k1)
    k3 = dt * f(t + 0.5*dt, y + 0.5*k2)
    k4 = dt * f(t + dt, y + k3)
    return t + dt, y + (k1 + 2*(k2 + k3) + k4)/6.0

inits = np.array([x for x in np.random.uniform(-1.5, 1.5, Nparticles*2)
    if x < -0.5 or x > 0.5][:Nparticles])

print(inits)

t0 = 0.0
tf = 2.0

plt.figure(figsize=(12,10))
plt.suptitle("Bohmian trajectories for %d particles" % len(inits), fontsize=20)
plt.xlabel("t")
plt.ylabel("y")

for index, y0 in enumerate(inits):
    
    x_list = []
    y_list = []
    y = y0
    t = t0
    
    while t < tf:
        t, yt = rk4(ydot, t, y, dt)
        x_list.append(kx*t)
        y_list.append(yt)
        y = yt
    
    plt.plot(x_list, y_list, color="black")
    plt.ylim(-13.0,13.0)
    filename = "TEMP_appendix_C-" + format("%05d" % index) + ".png"
    plt.savefig(filename, dpi=150)
    print("Saved to %s" % filename)
