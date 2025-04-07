#!/usr/bin/env python3
"""
Quantum mixed state evolver, and wigner distribution plotter.

Basic strategy:
- Use a spatial finite-difference scheme with not too many pixels. Just
  enough that pixelization doesn't become an issue (i.e., at least a few points
  per wavelength), but not so much that it becomes ridiculously slow.
- Just use full eigen-decomposition of the hamiltonian, don't bother with any
  time-stepping stuff.
- Create a density matrix by statistically superimposing various gaussian wave
  packet pure-states.
- Calculate density matrix in energy basis at any time in future, easy since
  the evolution operator is diagonal in this basis.
- Convert to wigner distribution by:
  - first transforming to the position basis,
  - then using a fourier transform trick (wigner distribution is just fourier
    transform of the antidiagonals!)
"""

from pylab import *

rc('text', usetex=True) ; rc('savefig', dpi=300)

### Parameters
# Very important number, smaller means more classical (finer-spaced discrete levels, larger means more quantum (fewer discrete levels)
hbar = 0.10/(2*pi)

def potential(x):
    return x**6 + 4*x**3 - 5*x**2 - 4*x
# the position-space:
x = linspace(-0.4,1.4,751)
N = len(x)
dx = x[1] - x[0]
U = potential(x)
mass = 1.0
# As an optimization we will throw away all high-energy eigenstates. They are going to
# be trash due to high spatial frequency (so finite-difference momentum is not
# accurate), and they just slow things down.
eig_cutoff = 100

### Calculate the Hamiltonian H = U + hbar^2 (d/dx)^2/2m
# Potential energy goes on diagonal of hamitonian
H = diag(U)
# For the kinetic energy operator (-hbar^2/(2*m) * d^2/dx^2), the usual central difference
# approximation will be used. These go on diagonal and off-diagonals.
H.flat[0::len(x)+1] += hbar*hbar/(mass*dx*dx) # on diagonal
H.flat[1::len(x)+1] += -0.5*hbar*hbar/(mass*dx*dx) # upper diagonal
H.flat[len(x)::len(x)+1] += -0.5*hbar*hbar/(mass*dx*dx) # lower diagonal

# Diagonalize the hamiltonian and cut off high energies
eigval, eigvec = eigh(H)
eigval = eigval[:eig_cutoff]
eigvec = eigvec[:,:eig_cutoff]
eigvec_H = conj(eigvec.T)
# This is basically the true hamiltonian we are going to model, with cutoff.
H_actual = eigvec @ diag(eigval) @ eigvec_H

### Prepare a density matrix (in energy eigenvector space) for time=0
rho_E = np.zeros((len(eigval), len(eigval)), dtype=complex)
# Use gaussian wavepacket pure states, then create a mixture with an uncertainty
# in the wavepacket position and momentum.
xwidth = sqrt(hbar)*0.33 ; countx = int(2 * 0.5 // xwidth)
pwidth = 0.5*hbar/xwidth ; county = int(2 * 2 // pwidth)
print(f"using {countx}*{county} wavepackets with sigma_x = {xwidth:.4f} and sigma_p = {pwidth:.4f}")
for x0 in linspace(0 + xwidth, 0.5 - xwidth, countx):
    for p0 in linspace(-1 + pwidth, 1 - pwidth, county):
        # make the packet and decompose it into eigenbasis
        packet = ((xwidth/dx)**2 * 2*pi)**(-0.25) * exp(-(x - x0)**2 / (2*xwidth)**2 + (1j * p0 / hbar)*x)
        decomp = eigvec_H @ packet
        # add it to matrix
        rho_E += decomp[:, None] * conj(decomp)
# normalize it
rho_E /= trace(rho_E)
assert np.all(rho_E == rho_E.T.conj()), "must be hermitian"

# It is straightforward to time-evolve rho_E

def rho_E_time(rho_E, t):
    # unitarily evolve the rho_E -- easy since it's the energy basis so
    # the unitary operator U = exp(-iHt/hbar) is diagonal.
    # rho(t) = U(t) rho(0) U(t)^h    so this is elementwise multiply by a phase shift.
    Ediff = eigval[:,None] - eigval
    return exp((-1j * t/hbar) * Ediff) * rho_E
# It's also straightforward to calculate rho_x (the position-basis density matrix).
def rho_x_from_E(rho_E):
    return eigvec @ rho_E @ eigvec_H

# Define some parameters of the wigner transform
wig_Np = int(1.5*N) # increase the ordinary momentum resolution by padding
wig_dp = hbar*pi/(dx*wig_Np) # p resolution
wig_p0 = -(wig_Np//2)*wig_dp
wig_extent = (x[0] - 0.5*dx,
              x[-1] + 0.5*dx,
              wig_p0 - 0.5*wig_dp,
              wig_p0 + (wig_Np + 0.5)*wig_dp)

def wignerify(rho_x):
    """ compute wigner transform of space-basis density matrix """
    # We want to FFT along the antidiagonals: [A] , [aBc], [bCd], and [D]
    # To do so we'll 'rotate the matrix' anticlockwise by 45 degrees
    #         ----        ----
    #        |Axax|      | ab |
    #        |xBxb|  --> |ABCD|
    #        |cxCx|      | cd |
    #        |xdxD|      |    |
    #         ----        ----
    # (this puts the diagonals onto rows)
    # The dropped 'x' elements contain mostly redundant information. They
    # could be included but that requires extra computation due to how they straddle
    # across the diagonal.
    Nx = len(rho_x)
    rhoflat = rho_x.ravel('C')
    rho_rot = np.zeros((Nx,Nx), dtype=complex, order='F')
    for i,row in enumerate(rho_rot):
        diagnum = (Nx - 1) // 2 - i
        if diagnum >= 0:
            row[diagnum:Nx-diagnum] = rhoflat[2*diagnum:Nx*(Nx+1-2*diagnum):Nx+1]
        else:
            row[-diagnum:Nx+diagnum] = rhoflat[-2*diagnum*Nx::Nx+1]

    # Fourier transform along the columns.
    res = fft(rho_rot, n = wig_Np, axis=0)
    # deshift the Nx // 2
    res *= exp((2j*pi*(Nx//2)/wig_Np) *  arange(wig_Np)[:,None])
    # Since input rho_x was hermitian, the wigner transform has no imaginary part
    # aside from numerical errors.
    res = res.real
    # put 0 momentum in middle and then snip to size
    res = np.roll(res, wig_Np//2, axis=0)
    return (1./(pi * hbar)) * res

def wigplot(wigner):
    fig = figure(figsize=(4,3))
    ax = axes((0.06,0.07,0.87,0.86))
    ax.set_xlim(-0.20,1.38)
    ax.set_ylim(-3.4,3.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('position $x$')
    ax.set_ylabel('momentum $p$')

    ax.imshow(wigner,
              origin='lower', extent=wig_extent,
              aspect='auto', cmap = matplotlib.cm.RdBu,
              interpolation='bicubic',
              norm=matplotlib.colors.SymLogNorm(linthresh=0.2, linscale=0.5,
                                                vmin=-4, vmax=4),
              )

    # add a rectangle to represent area of hbar
    hbarx = sqrt(hbar) * 0.40
    hbary = hbar/hbarx
    ax.add_artist(Rectangle((-0.07,-3.3),
                  hbarx*sqrt(2*pi), hbary*sqrt(2*pi), facecolor=(0.6,0.85,0.6), edgecolor='k', linewidth=0.5, linestyle='--'))
    ax.text(-0.07-0.01, -3.3 + 0.5*hbary*sqrt(2*pi), r'$h = {}\quad$', ha='right', va='center')

    ax.add_artist(Rectangle((-0.07,-2.4),
                  hbarx, hbary, facecolor=(0.6,0.85,0.6), edgecolor='k', linewidth=0.5, linestyle='--'))
    ax.text(-0.07-0.01, -2.4 + 0.5*hbary, r'$\hbar = {}\quad$', ha='right', va='center')

    # draw marginal x- and p- probability distributions on edges of plot
    ax_Px = axes((0.06, 0.93, 0.87, 0.07))
    ax_Px.axis('off')
    ax_Px.set_xlim(*ax.get_xlim())
    ax_Px.set_ylim(0,5.4)
    Px = np.sum(wigner, axis=0) * wig_dp
    ax_Px.fill_between(x, 0, Px, facecolor='#1f77b4')

    ax_Pp = axes((0.93, 0.07, 0.07, 0.86))
    ax_Pp.axis('off')
    ax_Pp.set_ylim(*ax.get_ylim())
    ax_Pp.set_xlim(0,1.0)
    Py = np.sum(wigner, axis=1) * dx
    ax_Pp.fill_betweenx(wig_p0 + arange(len(Py))*wig_dp, 0, Py, facecolor='#1f77b4')

    return fig, ax

def entropy(rho):
    rho_ev = eigvalsh(rho)
    return -sum(rho_ev * nan_to_num(log(rho_ev)))

#%%
# plot some eigenstates
def plot_eig(n):
    rho = zeros_like(rho_E) ; rho[n][n] = 1.
    assert entropy(rho) == 0.
    fig, ax = wigplot(wignerify(rho_x_from_E(rho_E_time(rho, 0))))
    fig.savefig(f'wigner_eigenstate_{n}.png')
plot_eig(0)
plot_eig(5)
plot_eig(30)
plot_eig(40)

#%%
# plot the equilibrated distribution
rho_E_equilib = rho_E * eye(len(rho_E))
fig, ax = wigplot(wignerify(rho_x_from_E(rho_E_time(rho_E_equilib, 0))))
fig.savefig('wigner_equilibrated.png')
# plot the initial distribution
fig, ax = wigplot(wignerify(rho_x_from_E(rho_E_time(rho_E, 0))))
fig.savefig('wigner_initial.png')
print(f"Initial entropy = {entropy(rho_E)}")
print(f"Equilibrated entropy = {entropy(rho_E_equilib)}")

#%%
# plot the significant (99.9%) Schmidt-decomposed pure states at times 0,0.5,1
def schmidt(rho_E, label):
    rho_eigval, rho_eigvec = eigh(rho_E)
    cumprob = 0.
    for i,(prob, vec) in enumerate(zip(rho_eigval[::-1], rho_eigvec.T[::-1])):
        cumprob += prob
        vec_x = eigvec @ vec # pure wavefunction in x space
        rho_x = vec_x[:,None] @ conj(vec_x)[None,:]
        fig, ax = wigplot(wignerify(rho_x))
        fig.text(0.5, 0.90, f'schmidt[{i}] {label} prob={prob:.3f}', ha='center', va='top')
        # plot
        ax.plot(x,10*vec_x.real+2., color='lime', ls='solid', lw=0.3, ms=1., mew=0.2, marker='+')
        ax.plot(x,10*vec_x.imag+2., color='lime', ls='dashed', lw=0.3)
        fig.savefig(f'schmidt/schmidt[{i}] {label}.png')
        close(fig)
        if cumprob > 0.999:
            break
schmidt(rho_E_time(rho_E,0.), label='t=0.0')
schmidt(rho_E_time(rho_E,0.5), label='t=0.5')
schmidt(rho_E_time(rho_E,1.), label='t=1.0')

# make a movie

framerate = 24
speed = 1.0 # how many time units per second
dt = speed / framerate

for i in range(0,741):
    if i <= 500:
        text = None
        t = i * dt
    elif i <= 620:
        # jump forward to a late time
        text = 'one hour later'
        t = 3600*speed + i * dt
    elif i <= 740:
        # jump forward to a late time
        text = 'one year later'
        t = 365.2425*86400*speed + i * dt
    wigner = wignerify(rho_x_from_E(rho_E_time(rho_E, t)))

    print(f"frame {i: 4d} range {amin(wigner):.03f} to {amax(wigner):.03f}")
    sys.stdout.flush()
    fig, ax = wigplot(wigner)
    if text:
        fig.text(0.5, 0.90, text, ha='center', va='top')
    fig.savefig(f'wigframes/frame{i:04d}.png')
    close(fig)

# commands to make the video file:
# cmd='ffmpeg -framerate 24 -f image2pipe -i - -crf 20 -b:v 0 -an -c:v libvpx-vp9'
# cat wigframes/*.png | $cmd -pass 1 -f webm -y /dev/null
# cat wigframes/*.png | $cmd -pass 2 -y out-full-20-vp9.webm