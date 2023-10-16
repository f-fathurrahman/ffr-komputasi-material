# modified by efefer

import numpy as np
# seed of random numbers is the system time by default

NITER = 100_000
n_in = 0
#################
### Main Loop ###
#################
for iiter in range(NITER):
    x = np.random.uniform(0, 1)
    y = np.random.uniform(0, 1)
    if x*x + y*y < 1.0:
        #n_in = n_in + 1
        n_in += 1
    #print(iiter+1, n_in/(iiter+1))

print("pi approx = ", n_in/(NITER)*4)


