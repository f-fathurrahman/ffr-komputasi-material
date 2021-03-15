import numpy as np
import matplotlib.pyplot as plt

ε0 = 2.0
t = 1.0
a = 1.0

# 1st BZ
k = np.linspace(-np.pi/a, np.pi/a, 100)
E = ε0 - 2*t*np.cos(k*a)

# Free electron
Efree = 0.5*k**2

plt.clf()
plt.plot(k, E, label="tight binding")
plt.plot(k, Efree, label="free electron")
plt.grid(True)
plt.legend()
plt.savefig("IMG_tb1d_monoatomik.pdf")
