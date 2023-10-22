import numpy as np
from numba import jit

@jit(nopython=True)
def do_sim(NITER):
    n_in = 0
    for i in range(NITER):
        x = np.random.uniform(0, 1)
        y = np.random.uniform(0, 1)
        if x*x + y*y < 1.0:
            n_in += 1
    pi_approx = 4 * n_in/NITER
    return pi_approx

# warm up
_ = do_sim(1)

import time
start = time.perf_counter()
pi_approx = do_sim(1_000_000)
end = time.perf_counter()
print("Elapsed (Numba) = {}s".format((end - start)))

print("pi_approx = %.10f" % pi_approx)