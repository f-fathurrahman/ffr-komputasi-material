Quantum ESPRESSO version used in this tutorial: 6.4

Calculating band structure is one of typical task in electronic structure calculations.
In this tutorial, I will describe how to calculate band structure of silicon crystal
using density functional theory as implemented in PWSCF which is included in the
Quantum EPSRESSO package.

We will use the following steps.

(1) Running an SCF calculation to find the converged Kohn-Sham potential

(2) Running a band structure calculation. This is simply solving
    Schroedinger equation with the Kohn-Sham potential obtained from
    the previous SCF calculation.

(3) Extract the band structure data and plot the result. 
    In this tutorial, I will use Python (with Matplotlib)
    to create the plot.

Let's start with the SCF calculation.
Here is the input file.

```
&CONTROL
  calculation = 'scf'
  pseudo_dir = '/home/efefer/pseudo'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 2
  celldm(1) = 10.20
  nat = 2
  ntyp = 1
  ecutwfc = 18.0
/

&ELECTRONS
  mixing_beta = 0.7
  conv_thr =  1.0d-8
/

ATOMIC_SPECIES
Si  28.086  Si.pz-vbc.UPF

ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25

K_POINTS automatic
  8 8 8 0 0 0
```

Let's name this file as `PWINPUT_scf` (you can use another name, of course).

Run this with the following command.

```
pw.x < PWINPUT_scf | tee LOG_scf
```

This calculation should not take too much time.

Next, we will do the band calculation (in the same directory as the SCF calculation).

Here is the input file. Let's name it `PWINPUT_bands`.

```
&CONTROL
  calculation = 'bands'
  pseudo_dir = '/home/efefer/pseudo'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 2
  celldm(1) = 10.20
  nat = 2
  ntyp = 1
  ecutwfc = 18.0
  nbnd = 8
/

&ELECTRONS
  mixing_beta = 0.7
  conv_thr =  1.0d-8
/

ATOMIC_SPECIES
Si  28.086  Si.pz-vbc.UPF

ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25

K_POINTS crystal
 60
0.50000000 0.25000000 0.75000000 0.00000000
0.50000000 0.28125000 0.71875000 0.15594042
0.50000000 0.31250000 0.68750000 0.31188083
0.50000000 0.34375000 0.65625000 0.46782125
0.50000000 0.37500000 0.62500000 0.62376166
0.50000000 0.40625000 0.59375000 0.77970208
0.50000000 0.43750000 0.56250000 0.93564249
0.50000000 0.46875000 0.53125000 1.09158291
0.50000000 0.50000000 0.50000000 1.24752332
0.47727273 0.47727273 0.47727273 1.38050976
0.45454545 0.45454545 0.45454545 1.51349619
0.43181818 0.43181818 0.43181818 1.64648262
0.40909091 0.40909091 0.40909091 1.77946905
0.38636364 0.38636364 0.38636364 1.91245549
0.36363636 0.36363636 0.36363636 2.04544192
0.34090909 0.34090909 0.34090909 2.17842835
0.31818182 0.31818182 0.31818182 2.31141479
0.29545455 0.29545455 0.29545455 2.44440122
0.27272727 0.27272727 0.27272727 2.57738765
0.25000000 0.25000000 0.25000000 2.71037409
0.22727273 0.22727273 0.22727273 2.84336052
0.20454545 0.20454545 0.20454545 2.97634695
0.18181818 0.18181818 0.18181818 3.10933339
0.15909091 0.15909091 0.15909091 3.24231982
0.13636364 0.13636364 0.13636364 3.37530625
0.11363636 0.11363636 0.11363636 3.50829268
0.09090909 0.09090909 0.09090909 3.64127912
0.06818182 0.06818182 0.06818182 3.77426555
0.04545455 0.04545455 0.04545455 3.90725198
0.02272727 0.02272727 0.02272727 4.04023842
0.00000000 0.00000000 0.00000000 4.17322485
0.02631579 0.00000000 0.02631579 4.30454309
0.05263158 0.00000000 0.05263158 4.43586134
0.07894737 0.00000000 0.07894737 4.56717958
0.10526316 0.00000000 0.10526316 4.69849783
0.13157895 0.00000000 0.13157895 4.82981607
0.15789474 0.00000000 0.15789474 4.96113432
0.18421053 0.00000000 0.18421053 5.09245256
0.21052632 0.00000000 0.21052632 5.22377081
0.23684211 0.00000000 0.23684211 5.35508905
0.26315789 0.00000000 0.26315789 5.48640729
0.28947368 0.00000000 0.28947368 5.61772554
0.31578947 0.00000000 0.31578947 5.74904378
0.34210526 0.00000000 0.34210526 5.88036203
0.36842105 0.00000000 0.36842105 6.01168027
0.39473684 0.00000000 0.39473684 6.14299852
0.42105263 0.00000000 0.42105263 6.27431676
0.44736842 0.00000000 0.44736842 6.40563501
0.47368421 0.00000000 0.47368421 6.53695325
0.50000000 0.00000000 0.50000000 6.66827150
0.50000000 0.04166667 0.54166667 6.81529353
0.50000000 0.08333333 0.58333333 6.96231556
0.50000000 0.12500000 0.62500000 7.10933760
0.50000000 0.16666667 0.66666667 7.25635963
0.50000000 0.20833333 0.70833333 7.40338166
0.50000000 0.25000000 0.75000000 7.55040370
0.46875000 0.28125000 0.75000000 7.66067022
0.43750000 0.31250000 0.75000000 7.77093675
0.40625000 0.34375000 0.75000000 7.88120327
0.37500000 0.37500000 0.75000000 7.99146980
```

There are several notable difference with the input file
as compared to the SCF input. The first one is that we have
used `calculation = 'bands'`.
The next, and probably the most trickiest part, is to choose the k-point
path. Usually, a band structure of a solid is plotted along the k-point
path. This k-point path is made by connecting several high-symmetry points
in the first Brilliouin zone. In the first example, I want to choose the k-point
path which is similar to the one used in this Figure:
[http://www.bandstructure.jp/Table/BAND/band_png/si_lda_5125.ps.png].

In several DFT package, we can simply define k-point path by specifying the
high-symmetry points and number of points to be sampled along each these
points. However (as far as I know), PWSCF does not have this feature, so,
we have to generate our k-point path by ourself. Fortunately, we can use
ASE to help us to generate this k-point path. An example of the script
I used to generate this k-point path is listed below.
This script will generate k-points path along
W-L-$\Gamma$-X-W-K.

```
from ase.dft.kpoints import *
from ase.units import Bohr
import sys

# Any lattice parameter should work, the important ones are the
# lattice vectors v1, v2, and v3. In this case we used the definition
# of FCC lattice vectors used in PWSCF.
alat = 10.20*Bohr
v1 = [-1,0,1]
v2 = [ 0,1,1]
v3 = [-1,1,0]
cell = 0.5*alat*np.transpose( np.array( [v1, v2, v3] ) )


# Number of total k-points in the path
NKPT = 60
# Use ase.dft module for obtaining k-points along high symmetry directions
points = ibz_points["fcc"]
W = points["W"]
L = points["L"]
G = points["Gamma"]
X = points["X"]
K = points["K"]
kpts, x, Xkpt = get_bandpath([W, L, G, X, W, K], cell, npoints=NKPT)

# Write kpts in the format understood by PWSCF
# The weights of the k-points are not used, so they can take any value.
# In this case we set them all to x[ik]
sys.stdout.write("K_POINTS crystal\n")
sys.stdout.write("%d\n" % NKPT)
for ik in range(NKPT):
    sys.stdout.write('%.8f %.8f %.8f %.8f\n' % (kpts[ik,0],kpts[ik,1],kpts[ik,2],x[ik]))
```

Now, we can run PWSCF code to calculate the band structure.
```
pw.x < PWINPUT_bands | tee LOG_bands
```

In the log file, you should see something similar to this:
```
     End of band structure calculation

          k =-1.0000 0.5000 0.0000 (   360 PWs)   bands (ev):

    -1.4205  -1.4205   2.2859   2.2859  10.4873  10.4873  11.2922  11.2922

          k =-0.9375 0.5000 0.0625 (   348 PWs)   bands (ev):

    -1.8314  -1.0316   2.0471   2.6556  10.1117  10.4752  11.3474  11.7794

          k =-0.8750 0.5000 0.1250 (   342 PWs)   bands (ev):

    -2.2236  -0.7384   2.0155   3.0987   9.7584  10.0893  11.9105  12.3847

    ... // snipped
```

This is the data for the band structure that we want. You can use your
favourite scripting language (or tools) to parse the log file and
build the data needed to be plotted as band structure.
However, if you don't want to do that we can use a utility program called
`bands.x` (also included in the Quantum ESPRESSO) to collect the band structure
data and produce an output which can be directly plotted. Quantum ESPRESSO
also includes an utility program to plot the band structure, however, in this
tutorial I will not use that and used Python to do the plotting.

Before plotting the band structure, let's see how to use `bands.x` to collect
the band structure data. We need an input file. Let's name this file as
`bands.inp`. The content of this file is:
```
&BANDS
  outdir = './tmp'
/
```
You should adjust this to your previous input files (SCF and bands calculation).

Run `bands.x` with the following command:
```
bands.x < bands.inp | tee LOG_pp_bands
```
There are various information (most related to symmetry) in the file `LOG_pp_bands`.
You will notice that after we run `bands.x`, several new files will be produced.
The file that we want to use is named `bands.out.gnu`.
This file consists of two columns. The first column is the coordinate
of k-points which is 'linearized'. The k-point coordinates in 3d dimension
is mapped into distances in 1d. We don't need to fuss about the details about this
in this case because it is already calculated by `bands.x` for us. 
The second column contains the band energies in eV. The data between
different bands are separated by one space.

In the file `LOG_pp_bands` we need the coordinates of special k-points.
This information is given in the following lines:
```
     high-symmetry point: -1.0000 0.5000 0.0000   x coordinate   0.0000
     high-symmetry point: -0.5000 0.5000 0.5000   x coordinate   0.7071
     high-symmetry point:  0.0000 0.0000 0.0000   x coordinate   1.5731
     high-symmetry point: -1.0000 0.0000 0.0000   x coordinate   2.5731
     high-symmetry point: -1.0000 0.5000 0.0000   x coordinate   3.0731
     high-symmetry point: -0.7500 0.7500 0.0000   x coordinate   3.4267
```
This is the coordinate of high-symmetry points that we specified when
we are generating the k-points path, i.e.:
W-L-G-X-W-K.


The script below can be used to plot the band structure data in the file
`bands.out.gnu`.

```
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
```

This the the result:
<figure>
  <img src="http://cmd.tf.itb.ac.id/wp-content/uploads/2019/07/Si_bands_v1.png"
       width="100%" height="100%">
</figure>

As exercise you may want to repeat the procedure for different k-points path
like in the following figures.

<figure>
  <img src="http://cmd.tf.itb.ac.id/wp-content/uploads/2019/07/Si_bands_v2.png"
       width="100%" height="100%">
</figure>


<figure>
  <img src="http://cmd.tf.itb.ac.id/wp-content/uploads/2019/07/Si_bands_v3.png"
       width="100%" height="100%">
</figure>


For more information about k-points path you can read the following:

[https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#special-points-in-the-brillouin-zone)].
