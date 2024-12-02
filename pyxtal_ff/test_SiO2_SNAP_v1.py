import time
from descriptors_SNAP import SO4_Bispectrum

from utilities_database import parse_OUTCAR_comp

data = parse_OUTCAR_comp("OUTCAR_comp")

lmax = 1
rcut = 4.0
w = {'Si': 2.0, 'O': 2.0}

start1 = time.time()
f = SO4_Bispectrum(w, lmax=lmax, rcut=rcut, derivative=False, stress=False, normalize_U=False, cutoff_function='cosine', rfac0=0.99363)
x = f.calculate(data[0]["structure"])
start2 = time.time()
print('time elapsed 1st: {}'.format(start2 - start1))

start1 = time.time()
f = SO4_Bispectrum(w, lmax=lmax, rcut=rcut, derivative=False, stress=False, normalize_U=False, cutoff_function='cosine', rfac0=0.99363)
x = f.calculate(data[0]["structure"])
start2 = time.time()
print('time elapsed 2nd: {}'.format(start2 - start1))
