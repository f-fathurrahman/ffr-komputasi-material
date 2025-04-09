import time

from my_pyxtal_ff.descriptors_SO3 import SO3
from my_pyxtal_ff.utilities_database import parse_OUTCAR_comp

data = parse_OUTCAR_comp("OUTCAR_comp")

lmax = 3
nmax = 4
rcut = 4.9
alpha = 2.0

start1 = time.time()
f = SO3(nmax=nmax, lmax=lmax, rcut=rcut, alpha=alpha, derivative=True, stress=False, cutoff_function='cosine')
x = f.calculate(data[0]["structure"])
start2 = time.time()
#print('x', x['x'])
#print('dxdr', x['dxdr'])
print('calculation time {}'.format(start2-start1))
