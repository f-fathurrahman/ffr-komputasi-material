import time
from math import pi
import numpy as np
import sympy
import matplotlib.pyplot as plt

def wavefunction(slit_distance, slit_width, st, velx, N, x, y, t):
    return N*sympy.exp((-((y-slit_distance)**2.0)/(4.0*slit_width*st)) + 1.0j*(velx*x-((t*velx**2.0)/2.0)))

def calc_psi(slit_distance, slit_width, st, velx, N, x, y, t):
    return wavefunction(slit_distance, slit_width, st, velx, N, x, y, t) + \
           wavefunction(slit_distance, slit_width, st, velx, N, x, -y, t)

# x is t
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

inits = np.array([x for x in np.random.uniform(-1.5, 1.5, Nparticles*2)
    if x < -0.5 or x > 0.5][:Nparticles])

print(inits)

t0 = 0.0
tf = 2.0

PSI = calc_psi(slit_distance, slit_width, st, velx, N, x, y, t)

sympy.pprint(PSI)
sympy.pprint(PSI.diff(y))

exit()

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
    filename = "TEMP_appendix_C-" + format("%05d" % index) + ".png"
    plt.savefig(filename, dpi=200)
    print("Saved to %s" % filename)
