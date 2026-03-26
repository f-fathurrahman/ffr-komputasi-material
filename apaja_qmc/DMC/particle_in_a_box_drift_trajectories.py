import numpy as np
import matplotlib.pyplot as plt
from numpy import pi, arccos, cos, exp
from matplotlib.patches import Rectangle
# https://matplotlib.org/stable/gallery/color/named_colors.html
from matplotlib.colors import CSS4_COLORS


rng = np.random.default_rng()
D = 1/2


# drift trajectory
def x_drift(tau, x0):
    return 1/pi*arccos(cos(pi*x0)*exp(-2*pi**2*D*tau))



# initial coordinates randomly in the unit box
Nw = 50 # number of walkers
dim = 2  # dimension
x0 = rng.uniform(size=(Nw,dim))

# plot trajectories
# =================

x = np.copy(x0)
plt.clf()
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.gca().set_aspect('equal') # square plot

Npoints = 100 # point on each trajectory
taus = np.linspace(0.0,0.2,Npoints) # stop trajectories at tau=0.2
xs = np.zeros((Nw,dim,Npoints))
for i,tau in enumerate(taus):
    xs[:,:,i] = x_drift(tau,x0)
    
    
# box:
box = Rectangle( (0.0,0.0), 1.0, 1.0, fill=None, color='blue', alpha=0.7 )
plt.gca().add_patch(box)
for i in range(Npoints):
    plt.scatter(xs[:,0,i],xs[:,1,i], c='maroon', s=0.5, alpha=0.4)
    plt.draw()
    plt.pause(1e-3)
    
# highlight start points as red
plt.scatter(xs[:,0,0],xs[:,1,0],c='red', s=1.0)
# highlight start points as blue
plt.scatter(xs[:,0,-1],xs[:,1,-1],c='blue', s=1.0)
plt.draw()
plt.pause(1e-3)

plt.show()





