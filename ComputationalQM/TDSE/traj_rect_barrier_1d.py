import numpy as np
import matplotlib.pyplot as plt

Nx = 1600  # Grid points

dt = 0.0001 # Evolution step

tmax = 6.0 # Propagation end

xmax = 50 # x- and y-window size
ymax = xmax

# 0 = periodic boundary
absorb_coeff = 20

# Initial wavefunction
# A Gaussian wave packet moving rightwards
def psi_0(x):
    vx = 10     # value of the initial velocity
    f = 0.j + np.exp(-((x+15)/4)**2)  # Gaussian profile
    f = f*np.exp(1.j*vx*x)   # Multiply by an x-dependent phase to introduce velocity
    return f

# A barrier modeled by V=0 for |x|>5 and V=40 for |x|<5
def eval_V(x, t, psi):
    V = np.piecewise(x, [abs(x-5)<2.5, abs(x-5)>=2.5],[40,0])
    return V;

def init_grid(Nx, xmax):
    x = np.linspace(-xmax, xmax-2*xmax/Nx, Nx)
    return x

# Builds the Laplacian in Fourier space
def build_Laplacian_fourier(Nx, xmax):
    kx = np.linspace(-Nx/4/xmax, Nx/4/xmax-1/2/xmax, Nx)     # x variable
    return (2*np.pi*1.j*kx)**2 

# Introduces an absorbing shell at the border of the computational window
def absorb_potential(x, xmax, dt, absorb_coeff):
    wx = xmax/20
    return np.exp(-absorb_coeff*(2-np.tanh((x+xmax)/wx)+np.tanh((x-xmax)/wx))*dt);

# Saves the data of abs(psi)**2 at different values of t
def save_psi2(psi):
    return abs(psi)**2



# builds spatial grid
x = init_grid(Nx, xmax)

# loads initial condition
psi = psi_0(x)

# Laplacian in Fourier space
L = build_Laplacian_fourier( Nx, xmax )

# linear phase in Fourier space (including point swap)
linear_phase = np.fft.fftshift( np.exp(1.j*L*dt/2) )

# Absorbing shell at the border of the computational window
border = absorb_potential(x, xmax, dt, absorb_coeff)

d_psi_real = np.gradient(np.real(psi), x)
d_psi_imag = np.gradient(np.imag(psi), x)

Nsteps = int(tmax/dt)

# potential operator
V = eval_V(x, 0, psi)

V_bohmian = np.imag( np.gradient(psi,x)/psi ) # no need to calculate velocity at this time?

from scipy.stats import rv_continuous

class PDF_psi0_gen(rv_continuous):

    def _pdf(self, x):
        #return np.exp(-(x+15)**2 / 2.) / np.sqrt(2.0 * np.pi)
        f = np.exp( -(x + 15.0)**2 /2.0 )/np.sqrt(2.0 * np.pi)
        return f #**2

# Initial positions
Nsamples = 10
initial_gen = PDF_psi0_gen(name="initial")
pos = initial_gen.rvs(size=Nsamples)

v = np.zeros(Nsamples)
dv = np.zeros(Nsamples)
dpos = np.zeros(Nsamples)

from scipy.interpolate import interp1d, splrep, splev, make_lsq_spline

# only include several points near xo
def interp_V_bohmian(xgrid, V, xo, dxgrid=0.5):
    x_min = xo - dxgrid
    x_max = xo + dxgrid
    idx = (xgrid >= x_min) & (xgrid <= x_max)
    #print("Using ", sum(idx), " for V_bohmian interpolations.")
    xx = xgrid[idx]
    f_interp = interp1d(xx, V[idx])
    return f_interp(xo)


# Main computational loop
for j in range(Nsteps):

    V[:] = eval_V(x, j*dt, psi)
    psi = psi*np.exp(-1.j*dt*V)
    psi = np.fft.fft(psi)
    psi = linear_phase*border*np.fft.ifft(psi)
    
    d_psi_real[:] = np.gradient(np.real(psi), x)
    d_psi_imag[:] = np.gradient(np.imag(psi), x)
    
    V_bohmian = np.imag( np.gradient(psi,x)/psi )
    
    # Negative velocity
    #f_interp = interp1d(x, V_bohmian)
    #v = f_interp( pos )

    # Not working
    #tck = splrep(x, V_bohmian)
    #v = splev(pos, tck, der=0)

    #t_spline = [-1, 0, 1]
    #k_spline = 3
    #t = np.r_[ (x[0],)*(k_spline+1), t_spline, (x[-1],)*(k_spline+1)]
    #f_interp = make_lsq_spline(x, V_bohmian, t_spline, k_spline)
    #v = f_interp( pos )

    for i in range(Nsamples):
        #v[i] = interp_V_bohmian(x, V_bohmian, pos[i])
        
        # RK4
        k1 = dt*interp_V_bohmian(x, V_bohmian, pos[i])
        k2 = dt*interp_V_bohmian(x, V_bohmian, pos[i] + 0.5*k1)
        k3 = dt*interp_V_bohmian(x, V_bohmian, pos[i] + 0.5*k2)
        k4 = dt*interp_V_bohmian(x, V_bohmian, pos[i] + k3)
        dpos[i] = (k1 + 2*(k2 + k3) + k4)/6.0


    #dpos = v*dt
    pos = pos + dpos
    t = (j+1)*dt
    
    #print("%18.10f %18.10f %18.10f %18.10f" % (t, v, dpos, pos))
    
    print("%18.10f" % t, end="")
    for i in range(Nsamples):
        print("; %18.10f" % pos[i], end="")
        #print("; %18.10f" % v[i], end="")
    print("")
