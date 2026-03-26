import numpy as np
import scipy as sc
from numpy import pi, cos, sin, exp, log, sqrt
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys

rng = np.random.default_rng()
dim = 3
D = 1/2
N = 2 # number of particles

E0_exact = N*dim*pi**2*D # 1 particle: 14.8044066, 2 bosons: 29.6088132,  etc.
ET = 10.0 # a guess
kappa = 10.0 # may need adjustment
  
Nwx = 20000
Nw_target = 5000
Nw = Nw_target

tau = 0.001 # default
if len(sys.argv)>1:
    tau = float(sys.argv[1])

plot = False # default
if len(sys.argv)>2:
    plot = int(sys.argv[2])

    
    
densfile =  'dens_2bosons_in_a_box_no_importance_sampling_tau='+str(tau)+'.dat'
print('saving density file ',densfile)

xw = np.zeros((Nwx,N,dim))
xw_tmp = np.zeros((Nwx,N,dim))       
alive = range(Nw)

# Energy error analysis
blocksize = int(1/(10*tau)) # just an approx
print('blocksize = ',blocksize)
E_data = 0.0
n_data = 0
nblocks = 100
E_bdata = np.zeros(nblocks)
n_bdata = 0

# walker update
def evolve(x):
    eta = rng.standard_normal(x.shape)
    d = np.sqrt(2*D*tau)*eta
    x += d
    return x

def branch(Nw, xw, alive):
    Mcopies = np.zeros(Nw,int)
    for i in alive:
        x = xw[i,:,:]
        if np.any(x<0) or np.any(x>1):
            # outside the box
            Mcopies[i] = 0.0
        else:
            xi = rng.uniform()
            Mcopies[i] = int(exp(tau*ET) + xi)

    ii = 0
    for i in alive:
        n = Mcopies[i]
        while n>0:            
            xw_tmp[ii] = np.copy(xw[i])
            ii += 1
            n -= 1
    Nw = ii
    alive = range(Nw)
    xw[alive] = xw_tmp[alive]
    return Nw, xw, alive
  


# Initial distribution
# ====================
xw[alive] = rng.uniform(low=0.1, high=0.9, size=(Nw,N,dim))

# The exact probability distribution function
# ===========================================
xs = np.linspace(0.0,1.0,1000)
# note that without importance sampling we get phi0, not phi0^2
# sin(pi*x) normalized to pdf
phi0 = pi/2*sin(pi*xs) * N

# Histogram parameters
# ====================
nbins = 100
bw = 1.0/nbins 
bins = np.linspace(0.0,1.0-bw,nbins) # leading edges

# histogram average counts
countscoll = np.zeros(nbins-1)
ncoll = 0


# DMC evolution loop
# ==================
t = 0.0
E0 = ET
nE = 1

measure = False


idmc = 0
while True:
    idmc += 1
    # DMC
    t += tau
    xw[alive] = evolve(xw[alive])
    Nw, xw, alive = branch(Nw, xw, alive)        
    ET = E0/nE + kappa*log(Nw_target/Nw)
    E0 += ET
    nE += 1
    if nE%10==0:
        print("t = %10.5f  <E0> = %15.8f  ET = %15.8f  Nw = %8d Exact E0 = %12.8f" %(t,E0/nE,ET,Nw,E0_exact))

    if measure:
        E_data += ET
        n_data += 1
        if n_data==blocksize:
            E_bdata[n_bdata] = E_data/n_data
            sigma = 0.0
            if n_bdata>1:
                sigma = np.std(E_bdata[:n_bdata])/sqrt(n_bdata-1)
            n_bdata += 1
            E_ave = E_bdata[:n_bdata].sum()/n_bdata
            print(80*'-')
            print('=== BLOCK ',n_bdata,'====')
            print("t = %10.5f  <E0> = %15.8f +/- %8.5f  Exact E0 = %12.8f" %(t,E_ave,sigma,E0_exact))
            print(80*'-')
            E_data = 0.0
            n_data = 0

            if n_bdata==nblocks:
                # end run
                break

        

    if plot:
        #plt.figure(3)
        #plt.title('Number of walkers')
        #ax = plt.gca()
        #ax.plot([nE],[Nw],'b.')
        #plt.draw()
        #plt.pause(1e-3)
        plt.figure(5)
        ax = plt.gca()
        ax.plot([t],[ET],'bo',markersize=1.)
        plt.draw()
        plt.pause(1e-3)
    
    # histogram of the probability distribution function

    if measure:
        x = xw[alive] 
        counts, _ = np.histogram(x, bins=bins) # combines all dimensions to x-axis
        countscoll += counts/(dim*Nw*bw)
        ncoll += 1
    
    if plot:
        x = xw[alive] 
        plt.figure(1)
        plt.clf()
        plt.xlim([-0.3,1.3])
        plt.ylim([-0.3,1.3])
        plt.gca().set_aspect('equal') # square plot
        plt.scatter(x[:,0],x[:,1],c='black', s=0.5)
        # box:
        box = Rectangle( (0.0,0.0), 1.0, 1.0, fill=None, color='blue', alpha=0.7 )
        plt.gca().add_patch(box)        
        plt.draw()
        plt.pause(1e-3)   
        
        plt.figure(2)
        plt.clf()
        plt.ylim([0,1.8])
        plt.title('instantaneous walker probability distribution')
        counts, _ = np.histogram(x, bins=bins) # combines all dimensions to x-axis
        plt.stairs(counts/(Nw*dim*bw), bins)
        plt.plot(xs, phi0, label= r'$\phi_0(x)$')
        plt.legend()    
        plt.draw()
        plt.pause(1e-3)

    # measurement 
    if t>0.2 and nE>100:
        if not measure:
            print(80*'=')
            print('Starting measurements')
            # restart energy measurements
            E0 = E0/nE            
            nE = 1
            measure = True

        if plot and ncoll>0:
            plt.figure(4)
            plt.clf()
            plt.ylim([0,1.8])
            plt.title('average walker probability distribution')
            plt.stairs(countscoll/ncoll, bins)
            plt.plot(xs, phi0, label= r'$\phi_0(x)$')
            plt.legend()
            plt.draw()
            plt.pause(1e-3)

    if measure and ncoll>0:
        dens = countscoll/ncoll
        with open(densfile,'w') as f:
            x = bw/2
            for d in dens:
                ex = pi/2*sin(pi*x)
                f.write('%10.8f  %10.8f  %10.8f \n' %(x,d,ex))
                x+=bw


Efile = 'E_2bosons_in_a_box_no_importance_sampling.dat'
print('adding result to file',Efile)
with open(Efile,'a') as f:
    f.write('%10.5f  %15.8f %15.8f %15.8f\n' %(tau,E0/nE,sigma,E0_exact))



if plot:
    plt.show()    
