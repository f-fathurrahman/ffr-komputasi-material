import numpy as np

# create x grid
dx = 0.1
x = np.arange(-10, 10+dx, dx)
N = len(x)

# create intial wave package
psi0 = np.exp(-0.25*x**2)/((2*np.pi)**0.25)

V = np.zeros([N])

# construct the 4th order FD matrix
g = -5j/(4*dx**2) - 1j*V
a = 1j/(24*dx**2)

diag_mat = np.diag(g)
off_diag1 = np.diag([16*a]*(N-1), 1) + np.diag([16*a]*(N-1), -1)
off_diag2 = np.diag([-a]*(N-2), 2) + np.diag([-a]*(N-2), -2)

M = diag_mat + off_diag1 + off_diag2

# time grid
dt = 0.01
t = np.arange(0, 20+dt, dt)
steps = len(t)

# create an array containing wavefunctions for each step
y = np.zeros([steps, N], dtype=np.complex128)
y[0] = psi0

import matplotlib.pyplot as plt
plt.clf()
plt.plot(x, np.abs(y[0]))
plt.savefig("TEMP_time_0.png", dpi=150)

# the RK4 method
for i in range(0, steps-1):
    k1 = np.dot(M, y[i])
    k2 = np.dot(M, y[i] + k1*dt/2)
    k3 = np.dot(M, y[i] + k2*dt/2)
    k4 = np.dot(M, y[i] + k3*dt)
    y[i+1] = y[i] + dt*(k1 + 2*k2 + 2*k3 + k4)/6

    plt.clf()
    plt.plot(x, np.abs(y[i+1]))
    plt.savefig("TEMP_time_" + str(i+1) + ".png", dpi=150)
