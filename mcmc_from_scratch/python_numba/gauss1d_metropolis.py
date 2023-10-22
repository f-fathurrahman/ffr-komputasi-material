from numba import jit
import numpy as np
from matplotlib import pyplot as plt
# seed of random numbers is the system time by default


@jit(nopython=True)
def do_sim(NITER):
    STEP_SIZE = 0.5

    # Set the initial configuration
    x = 0.0
    naccept = 0

    data_for_plot = []
    for i in range(NITER):
        backup_x = x
        action_init = 0.5*x*x
        # Random step size
        dx = np.random.uniform(-STEP_SIZE, STEP_SIZE)
        # Update
        x += dx
        action_fin = 0.5*x*x
        # Metropolis test
        metropolis = np.random.uniform(0,1)
        if(np.exp(action_init - action_fin) > metropolis):
            naccept += 1
        else:
            x = backup_x
        #print(x, naccept/(i+1))
        #
        data_for_plot.append(x)
    return data_for_plot


# Warm-up
_ = do_sim(1)


import time
start = time.perf_counter()
data_for_plot = do_sim(1_000_000)
end = time.perf_counter()
print("Elapsed (after compilation, Numba) = {}s".format((end - start)))


###########################################################
### plot the histogram and compare with analytic result ###
###########################################################
def gauss_func(x):
    return 1.0/np.sqrt(2.0*np.pi)*np.exp(-x*x/2.0)

y = np.arange(-4.0,4.0,0.1)
plt.xlabel("$x$")
plt.ylabel("$P(x)$")
plt.plot(y, gauss_func(y), label="exact")
plt.hist(data_for_plot, bins=80, density=True, label="MCMC")
plt.legend()
plt.show()


