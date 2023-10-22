import numpy as np
import os
import sys
from matplotlib import pyplot as plt
from numba import jit

# Calculation of the action
@jit(nopython=True)
def calc_action(spin, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY):
    action = 0.0
    sum1 = np.sum(spin)
    sum2 = 0
    for ix in range(NX):
        ixp1 = (ix+1)%NX #ixp1=ix+1; be careful about the boundary condition.
        for iy in range(NY):
            iyp1 = (iy+1)%NY #iyp1=iy+1; be careful about the boundary condition.
            #sum1=sum1+spin[ix,iy]
            sum2 += spin[ix,iy]*spin[ixp1,iy]+spin[ix,iy]*spin[ix,iyp1]
    action = (sum2*COUPLING_J + sum1*COUPLING_h)/TEMPERATURE*(-1.0)
    return action


# Calculation of the change of the action when the spin at (ix,iy) is flipped
@jit(nopython=True)
def calc_action_change(spin, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY, ix, iy):
    action_change = 0.0

    ixp1 = (ix+1)%NX #ixp1=ix+1; be careful about the boundary condition.
    iyp1 = (iy+1)%NY #iyp1=iy+1; be careful about the boundary condition.
    ixm1 = (ix-1+NX)%NX #ixm1=ix-1; be careful about the boundary condition.
    iym1 = (iy-1+NY)%NY #iym1=iy-1; be careful about the boundary condition.
    #
    sum1_change = 2*spin[ix,iy]
    sum2_change = 2*spin[ix,iy]*spin[ixp1,iy] + \
                  2*spin[ix,iy]*spin[ix,iyp1] + \
                  2*spin[ix,iy]*spin[ixm1,iy] + \
                  2*spin[ix,iy]*spin[ix,iym1]
    # 
    action_change = (sum2_change*COUPLING_J + sum1_change*COUPLING_h)/TEMPERATURE
    #
    return action_change

# Calculation of the total spin
@jit(nopython=True)
def calc_total_spin(spin):
    return np.sum(spin)


def init_config(NCONFIG):
    # Set the initial configuration
    if NCONFIG == 1:
        spin = np.ones((NX,NY))
    elif NCONFIG == -1:
        spin = np.ones((NX,NY))
        spin = spin*(-1)
    elif NCONFIG == 0:
        if os.path.exists('input_config.txt'):
            spin = np.empty((NX,NY))
            for l in open('input_config.txt').readlines():
                read = l[:-2].split(' ')
                ix = int(l[0])
                iy = int(l[1])
                spin[ix,iy] = read[2]
        else:
            print('no input configuration')
            sys.exit()
    return spin



@jit(nopython=True)
def do_sim(spin, NITER, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY):
    naccept = 0 #counter for the number of acceptance
    NSKIP = 1000
    for i in range(NITER):
        #choose a point randomly.
        rand_site = np.random.uniform(0, NX*NY)
        ix = int(rand_site/NX)
        iy = int(rand_site%NY) # huh?
        action_change = calc_action_change(spin, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY, ix, iy)
        #
        metropolis = np.random.uniform(0,1)
        if np.exp(-action_change) > metropolis:
            # accept
            spin[ix,iy] = -spin[ix,iy]
            naccept += 1
        # Reject: no change

        total_spin = calc_total_spin(spin)
        energy = calc_action(spin, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY)*TEMPERATURE

        if (i+1)%NSKIP == 0:
            print(f"Iteration = {i+1} is done")
        
    return spin




NITER = 4_096_000
NX = 64 #number of sites along x-direction
NY = 64 #number of sites along y-direction
COUPLING_J = 1.0
COUPLING_h = 0.1
TEMPERATURE = 5.0

NCONFIG = 1 #0 -> read 'input_config.txt'; 1 -> all up; -1 -> all down

# Warm up
spin = init_config(NCONFIG)
_ = do_sim(spin, 1, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY)


spin = init_config(NCONFIG)
import time
start = time.perf_counter()
spin = do_sim(spin, NITER, COUPLING_J, COUPLING_h, TEMPERATURE, NX, NY)
end = time.perf_counter()
print("Elapsed (Numba) = {}s".format((end - start)))


"""
######################################
### 2D plot of final configuration ###
######################################
plt.figure()
plt.imshow(spin, interpolation='nearest', vmin=0, vmax=1)
plt.show()
            
#########################
### save final config ###
#########################
path_output_config = 'TEMP_output_config.txt'
with open(path_output_config, mode='w') as output_config:
    for ix in range(NX):
        for iy in range(NY):
            print(ix,iy,spin[ix,iy], file=output_config)
"""
