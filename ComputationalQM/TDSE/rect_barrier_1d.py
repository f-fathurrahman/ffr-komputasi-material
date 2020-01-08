# A Gaussian wave packet impinges on a rectangular barrier. Part of the wave
# is reflected, part of the wave is transmitted.

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.style
matplotlib.style.use("dark_background")

# Grid points
Nx = 1600 

# Evolution step
dt = 0.0001

# Propagation end
tmax = 6

# x- and y-window size
xmax = 50
ymax = xmax

# number of .png images
images = 100

# 0 = periodic boundary
absorb_coeff = 20

# If 1, it plots on the screen but does not save the images
# If 2, it saves the images but does not plot on the screen
# If 3, it saves the images and plots on the screen
output_choice = 3

# Fixes a maximum scale of |psi|**2 for the plots. If 0, it does not fix it.
fixmaximum = 1.05

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

# Number of computational steps between consecutive graphic outputs
steps_image = int(tmax/dt/images)
Nsteps = steps_image*images + 1
idx_image = 1


# Main computational loop
for j in range(Nsteps):        # propagation loop
    
    print("step = ", j)
    
    # Generates image output
    if j % steps_image == 0:
        #plt.clf()
        #plt.plot(x, np.abs(psi)**2, label="Psi0^2")
        #plt.plot(x, border, label="border")
        #plt.ylim([0.0, 2.0])
        #plt.grid()
        #plt.legend()
        #plt.savefig("IMG_amp_psi0_{:08d}.png".format(idx_image), dpi=150)
        #
        plt.clf()
        plt.plot(x, d_psi_real, label="Re")
        plt.plot(x, d_psi_imag, label="Im")
        plt.ylim([-10.0, 10.0])
        plt.grid()
        plt.legend()
        plt.savefig("IMG_d_psi_{:08d}.png".format(idx_image), dpi=150)
        #
        idx_image = idx_image + 1


    V = eval_V(x, j*dt, psi)            # potential operator
    psi *= np.exp(-1.j*dt*V)            # potential phase
    psi = np.fft.fft(psi)            # 1D Fourier transform
    psi *=linear_phase                # linear phase from the Laplacian term
    psi = border*np.fft.ifft(psi)    # inverse Fourier transform and damping by the absorbing shell
    
    d_psi_real[:] = np.gradient(np.real(psi), x)
    d_psi_imag[:] = np.gradient(np.imag(psi), x)
