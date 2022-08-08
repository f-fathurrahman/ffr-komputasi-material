import random
from math import pi

def direct_pi(N):
    Nhits = 0
    for i in range(N):
        x = random.uniform(-1.0, 1.0)
        y = random.uniform(-1.0, 1.0)
        if x ** 2 + y ** 2 < 1.0:
            Nhits += 1
    return Nhits
 
Nruns = 1000
Ntrials = 1_000_000
for irun in range(Nruns):
    Nhits = direct_pi(Ntrials)
    pi_approx = 4.0 * Nhits / Ntrials
    error = abs(pi - pi_approx)
    print("%5d  %18.10f  %10.5e" % (irun, pi_approx, error))
