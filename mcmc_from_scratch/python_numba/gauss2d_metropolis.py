import numpy as np
from matplotlib import pyplot as plt

from numba import jit

@jit(nopython=True)
def do_sim(NITER):
    STEP_SIZE_x = 0.5
    STEP_SIZE_y = 0.5
    # Set the initial configuration ###
    x = 0.0
    y = 0.0
    naccept = 0

    # data for plot
    data_x = []
    data_y = []

    # Main Loop
    for i in range(NITER):
        backup_x = x
        backup_y = y
        action_init = 0.5*(x*x + y*y - x*y)
        dx = np.random.uniform(-STEP_SIZE_x, STEP_SIZE_x)
        dy = np.random.uniform(-STEP_SIZE_y, STEP_SIZE_y)
        #
        x += dx
        y += dy
        action_fin = 0.5*(x*x + y*y - x*y)

        # Metropolis test ###
        metropolis = np.random.uniform(0,1)
        if(np.exp(action_init-action_fin) > metropolis):
            naccept += 1
        else:
            # reject
            x = backup_x
            y = backup_y
    
        # data output
        #if (i+1)%10 == 0:
        #    #output the results every ten steps.   
        #    #print(x,y,naccept/(iter+1))
        
        data_x.append(x)
        data_y.append(y)
    #
    return data_x, data_y


# warm up
_, _ = do_sim(1)


import time
start = time.perf_counter()
data_x, data_y = do_sim(1_000_000)
end = time.perf_counter()
print("Elapsed (Numba) = {}s".format((end - start)))


# make 2D plot

# contour of the target distribution is plotted by white line
z1 = np.arange(-5, 5, 0.1)
z2 = np.arange(-5, 5, 0.1)
X,Y = np.meshgrid( z1, z2 )
Z = np.exp(-(X**2 + Y**2 - X*Y)/2) # analytical distrib
plt.contour(X, Y, Z,colors=['white'])

# plot 2D histogram
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.hist2d(data_x, data_y, bins=100, density=True)
plt.show()
