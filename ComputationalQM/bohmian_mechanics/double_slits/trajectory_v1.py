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

"""
y = sympy.var("y")
t = sympy.var("t")
x = y*t

st = sigma_0*( 1.0 + (1.0j*t)/(2*sigma_0**2) )
N = (2*pi*st**2)**-0.25

wavefunc_analytic = N*sympy.exp((-((y - Y)**2.0)/(4.0*sigma_0*st)) + 1.0j*(kx*x - ((t*kx**2.0)/2.0)))

print(wavefunc_analytic)

ydot_analytic = wavefunc_analytic.diff(y)/wavefunc_analytic
print(sympy.simplify(ydot_analytic))
"""


I = 1.0j

def ydot(t, y):
    yy = ( Y*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t) ) - Y + \
         2*I*sigma_0**2*t*v_x*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) + \
         2*I*sigma_0**2*t*v_x - t**2*v_x*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) - \
         t**2*v_x - y*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) - y)/(2*sigma_0**2*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) + \
         2*sigma_0**2 + I*t*cmath.exp(2*Y*y/(2*sigma_0**2 + I*t)) + I*t)
    #return np.imag(yy) + np.sign( np.random.randn() )
    return np.imag(yy) + 0.5*np.random.randn()

def xdot():
    return v_x

def rk4(x, h, y, f):
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5*h, y + 0.5*k1)
    k3 = h * f(x + 0.5*h, y + 0.5*k2)
    k4 = h * f(x + h, y + k3)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

def trajectories(t, state):
    y, vel = state
    return ydot(t, y)

#def rk4(f, t, y, dt):
#    k1 = dt * f(t, y)
#    k2 = dt * f(t + 0.5*dt, y + 0.5*k1)
#    k3 = dt * f(t + 0.5*dt, y + 0.5*k2)
#    k4 = dt * f(t + dt, y + k3)
#    return t + dt, y + (k1 + 2*(k2 + k3) + k4)/6.0

inits = np.array([x for x in np.random.uniform(-1.5, 1.5, Nparticles*2)
    if x < -0.5 or x > 0.5][:Nparticles])

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
