import numpy as np

def do_sim(NITER):
    sum_y = 0.0
    for iiter in range(NITER):
        x = np.random.uniform(0, 1)
        y = np.sqrt(1.0-x*x)
        sum_y += y
    pi_approx = 4 * sum_y/NITER
    return pi_approx

NITER = 1_000_000
import time
start = time.perf_counter()
pi_approx = do_sim(NITER)
end = time.perf_counter()
print("Elapsed (standard) = {}s".format((end - start)))

print("NITER = %d" % NITER)
print("pi approx = ", pi_approx)
print("diff = %e" % abs(pi_approx - np.pi))