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
points = ibz_points['fcc']
G = points['Gamma']
G1 = [1.0, 1.0, 1.0]
X = points['X']
L = points['L']
kpts, x, Xkpt = get_bandpath([L, G, X, G1], cell, npoints=NKPT)

# Write kpts in the format understood by PWSCF
# The weights of the k-points are not used, so they can take any value.
# In this case we set them all to 1.0
sys.stdout.write("K_POINTS crystal\n")
sys.stdout.write("%d\n" % NKPT)
for ik in range(NKPT):
  sys.stdout.write('%.8f %.8f %.8f %.8f\n' % (kpts[ik,0],kpts[ik,1],kpts[ik,2],x[ik]))

