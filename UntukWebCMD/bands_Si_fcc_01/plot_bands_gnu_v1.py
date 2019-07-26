import numpy as np
import matplotlib.pyplot as plt

NBANDS   = 8
NKPOINTS = 60

databands = np.loadtxt("bands.out.gnu_v1")

Nocc = 4  # number of occupied bands
PLOT_SAVE = "Si_bands_v1.pdf"
FIGSIZE = (6, 8)

# obtained from the output log of bands.x
Xkpt = [0.0000, 0.7071, 1.5731, 2.5731, 3.0731, 3.4267]
labelX = ["W", "L", r"$\Gamma$", "X", "W", "K"]

ebands = np.zeros( (NKPOINTS, NBANDS) )
kvec   = np.zeros( (NKPOINTS, NBANDS) )

for ib in range(NBANDS):
    idx1 = (ib)*NKPOINTS
    idx2 = (ib+1)*NKPOINTS
    ebands[:,ib] = databands[idx1:idx2,1]
    kvec[:,ib]   = databands[idx1:idx2,0]

Efermi = np.max( ebands[:,Nocc-1] )
ebands[:,:] = ebands[:,:] - Efermi

Emin = np.min(ebands)
Emax = np.max(ebands)

plt.figure(figsize=FIGSIZE)
plt.clf()
for ib in range(NBANDS):
    plt.plot( kvec[:,ib], ebands[:,ib] )

for p in Xkpt:
    plt.plot([p, p], [Emin, Emax], "k-")
plt.xticks(Xkpt, labelX)

plt.plot([0, Xkpt[-1]], [0, 0], "k--")
plt.ylim( Emin, Emax )
plt.xlim( 0, Xkpt[-1] )

plt.xlabel("k vector")
plt.ylabel("Energy (eV)")
plt.title("Band structure of Si")
plt.savefig(PLOT_SAVE)
plt.savefig("Si_bands_v1.png", dpi=200)

