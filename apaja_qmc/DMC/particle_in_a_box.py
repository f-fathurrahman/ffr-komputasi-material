import numpy as np
from numpy import pi, cos, sin, exp
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

rng = np.random.default_rng()
dim = 3
D = 1/2
Ntry = 0
Nacc = 0


# trial wf 
def varphi_T(x):
    if np.any(x<0) or np.any(x>1):
        res = 0.0
    else:
        res = np.prod(sin(pi*x)) 
    return res

# G(x' <- x; tau)
def G(xp, x):
    return exp(-np.sum((xp-x-D*tau*F(x))**2)/(4*D*tau))


# drift
def F(x):
    if np.any(x<0) or np.any(x>1):
        F = 0.0*x
    else:
        F = 2*pi*(cos(pi*x)/sin(pi*x))
    return F

# walker update
def evolve(x, Nacc, Ntry):
    eta = rng.standard_normal(x.shape)
    d = D*tau*F(x) + np.sqrt(2*D*tau)*eta
    xp = x + d

    # UNCOMMENT NEXT LINE TO SKIP THE ACCEPT/REJECT STEP
    #return xp, Nacc+1, Ntry+1 
    
    # Accept/reject step
    for i in range(Nw):
        xi = x[i,:]
        xpi = xp[i,:]
        Ntry += 1
        ratio = varphi_T(xpi)**2*G(xi,xpi)
        ratio/= varphi_T(xi)**2*G(xpi,xi)
        # Metropolis
        if ratio >= 1.0 :
            Nacc += 1
            x[i,:] = xpi
        else:
            r = rng.uniform()
            if ratio >= r:
                Nacc += 1
                x[i,:] = xpi
                
    return x, Nacc, Ntry


Nw = 5000
tau = 0.05 

# Initial distribution
# ====================
# choose one:
# start with
# 1)  walkers in the middle of the box
x = np.zeros((Nw,dim))+0.5
# 2) walkers scattered randomly in the box 
#x = rng.uniform(size=(Nw,dim))
# 3) walkers scattered randomly in the box avoiding edges 
#x = rng.uniform(low=0.1, high=0.9, size=(Nw,dim))


# The exact probability distribution function
# ===========================================
xs = np.linspace(0.0,1.0,1000)
phi0s = 2*sin(pi*xs)**2 # normalized pdf for (sin(pi*x))^2

#print(phi0s.sum()*(xs[1]-xs[0]))
#exit()


# Histogram parameters
# ====================
nbins = 101 
bins = np.linspace(0.0,1.0,nbins)
bw = bins[1]-bins[0]
# histogram average counts
countscoll = np.zeros(nbins-1)
ncoll = 0


# DMC evolution loop
# ==================
t = 0.0
while True:
    
    # walker positions are now possible positions of the particle
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


    # evolve for a while
    for j in range(10):
        t += tau
        x, Nacc, Ntry  = evolve(x, Nacc, Ntry)
        print(f't = {t:.7f}')
        
    # histogram of the probability distribution function
    plt.figure(2)
    plt.clf()
    plt.title('instantaneous walker probability distribution')
    counts, _ = np.histogram(x, bins) # combines all dimensions to x-axis
    plt.stairs(counts/(Nw*dim*bw), bins)    
    plt.plot(xs, (2**(1/2)*sin(pi*xs))**2, label= r'$\phi_0(x)^2$')
    plt.legend()    
    plt.draw()
    plt.pause(1e-3)


    # measurement 
    if t>0.1:
        countscoll += counts
        ncoll += 1
        plt.figure(3)
        plt.clf()
        plt.title('average walker probability distribution')
        plt.stairs(countscoll/(Nw*dim*bw*ncoll), bins)
        plt.plot(xs, (2**(1/2)*sin(pi*xs))**2, label= r'$\phi_0(x)^2$')
        plt.legend()
        plt.draw()
        plt.pause(1e-3)

    if Ntry>0:
        print(f'acceptance = {Nacc/Ntry*100.0:.6f}')


plt.show()    
