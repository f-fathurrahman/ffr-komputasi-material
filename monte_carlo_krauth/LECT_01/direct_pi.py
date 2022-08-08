import random
from math import pi

Ntrials = 1_000_000
Nhits = 0
for iiter in range(Ntrials):
    x = random.uniform(-1.0, 1.0)
    y = random.uniform(-1.0, 1.0)
    if x**2 + y**2 < 1.0: 
        Nhits += 1

pi_approx = 4.0 * Nhits / Ntrials
error = abs(pi - pi_approx)
print("pi_approx = %18.10f error = %10.5e" % (pi_approx, error))
