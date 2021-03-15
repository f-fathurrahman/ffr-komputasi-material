from math import sqrt, sin, cos

κ1 = 1.0
κ2 = 1.1
m = 1.0
a = 1.0

def real_wave(ω, k, t, n):
    return cos(ω*t - k*n*a)

def calc_omega(k):
    k1p2 = κ1 + κ2
    k1t2 = κ1*κ2
    ω1 = sqrt( (κ1 + κ2)/m + (1.0/m)*sqrt( (κ1 + κ2)**2 - 4*κ1*κ2*sin(k*a/2)**2 ) )
    ω2 = sqrt( (κ1 + κ2)/m - (1.0/m)*sqrt( (κ1 + κ2)**2 - 4*κ1*κ2*sin(k*a/2)**2 ) )
    return ω1, ω2

A1 = 1.0
A2 = 1.0
k = 2.0
ω1, ω2 = calc_omega(k)
print("oemga = ", ω1, ω2)

import matplotlib.pyplot as plt

# Equilibrium pos
X1 = -2.0
X2 =  2.0

Δt = 0.05
for i in range(100):
    t = i*Δt
    x1 = A1*real_wave(ω2, k, t, 1)
    x2 = A2*real_wave(ω2, k, t, 1)
    print("%18.10f %18.10f" % (x1, x2))
    plt.clf()
    plt.plot([X1+x1, X2+x2], [0.0, 0.0], marker="o", markersize=10, linewidth=0)
    plt.ylim([-1.0, 1.0])
    plt.xlim([-5.0, 5.0])
    filename = "IMG_t_%04d.png" % i
    plt.savefig(filename, dpi=100)

