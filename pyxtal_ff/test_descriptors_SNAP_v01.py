from ase.io import read
import time

from descriptors_SNAP import SO4_Bispectrum

from ase.build import bulk
test = bulk('Si', 'diamond', a=5.459)
cell = test.get_cell()
cell[0,1] += 0.5
test.set_cell(cell)

print(test)

lmax = 1
rcut = 4.0
w = {'Si': 2.0}

start1 = time.time()
f = SO4_Bispectrum(w, lmax=lmax, rcut=rcut, derivative=False, stress=False, normalize_U=False, cutoff_function='cosine', rfac0=0.99363)
x = f.calculate(test)
start2 = time.time()
print('time elapsed 1st: {}'.format(start2 - start1))

start1 = time.time()
f = SO4_Bispectrum(w, lmax=lmax, rcut=rcut, derivative=False, stress=False, normalize_U=False, cutoff_function='cosine', rfac0=0.99363)
x = f.calculate(test)
start2 = time.time()
print('time elapsed 2nd: {}'.format(start2 - start1))

#for key, item in x.items():
#    print(key, item)
#print(x['rdxdr'].shape)
#print(x['rdxdr'])
#print(np.einsum('ijklm->klm', x['rdxdr']))
print(x['x'])
