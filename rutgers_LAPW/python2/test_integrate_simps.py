from __future__ import print_function
import scipy.integrate
import numpy as np

x = np.linspace(0.0, 1.0, 201)
f = 3*np.cos(x)

s = scipy.integrate.simps(f, x)
print("s = ", s)
