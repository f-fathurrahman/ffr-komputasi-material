import time
from math import pi
import cmath
import numpy as np
import sympy
import matplotlib.pyplot as plt

Nparticles = 200
Y = 1.0
sigma_0 = 0.2
v_x = 0.1
dt = 1e-2

I = 1.0j

def ydot(t, y):
    yy = ( Y*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t) ) - Y + \
         2*I*sigma_0**2*t*v_x*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) + \
         2*I*sigma_0**2*t*v_x - t**2*v_x*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) - \
         t**2*v_x - y*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) - y)/(2*sigma_0**2*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) + \
         2*sigma_0**2 + I*t*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) + I*t)
    return np.imag(yy) + 10.0*np.random.randn()  # tambahkan Budiyono-Rohrlich term

def xdot():
    return v_x

def rk4(t, h, y, f):
    k1 = h * f(t, y)
    k2 = h * f(t + 0.5*h, y + 0.5*k1)
    k3 = h * f(t + 0.5*h, y + 0.5*k2)
    k4 = h * f(t + h, y + k3)
    return t + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

def trajectories(t, state):
    y, vel = state
    return ydot(t, y)

inits = np.array([x for x in np.random.uniform(-1.5, 1.5, Nparticles*2)
    if x < -0.5 or x > 0.5][:Nparticles])

# Formalnya: generate or sampling from P(x) = |psi(x)|^2


#print(inits)

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
    
    state = np.array([y0,0])
    while t < tf:
        t, state = rk4(t, dt, state, trajectories)
        x_list.append(t)
        y_list.append(state[0])

    #while t < tf:
    #    t, yt = rk4(ydot, t, y, dt)
    #    x_list.append(v_x*t)
    #    y_list.append(yt)
    #    y = yt
    
    plt.plot(x_list, y_list, color="black")
    plt.ylim(-13.0,13.0)
    filename = "TEMP_appendix_C-" + format("%05d" % index) + ".png"
    plt.savefig(filename, dpi=150)
    print("Saved to %s" % filename)
