import time
from math import pi, sqrt, sin, cos
import numpy as np
import sympy
import matplotlib.pyplot as plt

def wavefunction(y, t):
    # initial conditions
    a = 3.
    c = 1.
    k0 = sqrt(2.)*pi/2.
    x = 1.
    sint = sin(pi/4.)
    cost = cos(pi/4.)
    return ( sympy.exp( 1.j*(k0*(x*cost-y*sint)-c*t) ) *
             sympy.exp( -((x*cost-y*sint)-c*t)**2/4/a**2 ) *
             sympy.exp( -(x*sint+y*cost)**2/4/a**2 ) +
             sympy.exp( 1.j*(k0*(x*cost+y*sint)-c*t) ) *
             sympy.exp( -((x*cost+y*sint)-c*t)**2/4/a**2 ) *
             sympy.exp( -(x*sint-y*cost)**2/4/a**2) )

def calc_psi(y, t):
    return wavefunction(y, t)

def rk4(wave, del_wave, x, h, y, f):
    k1 = h * f(wave, del_wave, x, y)
    k2 = h * f(wave, del_wave, x + 0.5*h, y + 0.5*k1)
    k3 = h * f(wave, del_wave, x + 0.5*h, y + 0.5*k2)
    k4 = h * f(wave, del_wave, x + h, y + k3)
    return x + h, y + (k1 + 2*(k2 + k3) + k4)/6.0

def trajectories(wave, del_wave, t, state):
    y, vel = state
    return np.imag(del_wave(y,t)/wave(y,t))



start = time.time()

Nparticles = 200
slit_distance = 1.0
slit_width = 0.2
velx = 0.1
dt = 1e-2

y = sympy.var("y")
t = sympy.var("t")
x = y*t

st = slit_width*( 1 + (1j*t)/(2*slit_width**2) )
N = (2*pi*st**2)**-0.25

inits = np.linspace(-10.0, -30.0, 50)
inits = np.append(inits, np.linspace(10.5, 30.5, 50) )

print(inits)

t0 = -15.0
tf = 15.0

PSI = calc_psi(y, t)

print(PSI)

wave = sympy.lambdify([y,t], PSI,"numpy")
del_wave = sympy.lambdify([y,t], PSI.diff(y),"numpy")


plt.figure(figsize=(12,10))
plt.xlim(t0, tf)
plt.suptitle("Bohmian trajectories for %d particles" % len(inits), fontsize=20)
plt.xlabel("t")
plt.ylabel("y")

for index,y0 in enumerate(inits):
    x = []
    y = []
    t = t0
    state = np.array([y0,0])
    while t < tf:
        t, state = rk4(wave, del_wave, t, dt, state, trajectories)
        x.append(t)
        y.append(state[0])
    plt.plot(x,y,color="black")
    plt.ylim(-13.0,13.0)
    print("Completed trajectory %d completed" % index)
    print("Elapsed time: %s seconds" % str(time.time() - start))
    filename = "TEMP_appendix_F_-" + format("%05d" % index) + ".png"
    plt.savefig(filename, dpi=150)
    print("Saved to %s" % filename)
