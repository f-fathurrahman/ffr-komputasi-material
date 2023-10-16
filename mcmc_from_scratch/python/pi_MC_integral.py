import numpy as np
# seed of random numbers is the system time by default

NITER = 1_000_000
#NITER = 100_000
sum_y = 0.0
for iiter in range(NITER):
    x = np.random.uniform(0, 1)
    y = np.sqrt(1.0-x*x)
    sum_y += y
    #print(iiter+1, sum_y/(iiter+1))

pi_approx = sum_y/(NITER)*4
print("NITER = %d" % NITER)
print("pi approx = ", pi_approx)
print("diff = %e" % abs(pi_approx - np.pi))