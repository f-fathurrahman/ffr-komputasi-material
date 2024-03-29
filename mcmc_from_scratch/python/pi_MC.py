# modified by efefer

import numpy as np
# seed of random numbers is the system time by default

def do_sim(NITER):
    n_in = 0
    for i in range(NITER):
        x = np.random.uniform(0, 1)
        y = np.random.uniform(0, 1)
        if x*x + y*y < 1.0:
            n_in += 1
    pi_approx = 4 * n_in/NITER
    return pi_approx


import time
start = time.perf_counter()
pi_approx = do_sim(1_000_000)
end = time.perf_counter()
print("Elapsed (standard) = {}s".format((end - start)))

print("pi_approx = %.10f" % pi_approx)