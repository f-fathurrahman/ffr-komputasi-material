from ase.io import read
import time
from optparse import OptionParser

from descriptors_SO4 import SO4_Bispectrum

# ---------------------- Options ------------------------
parser = OptionParser()
parser.add_option("-c", "--crystal", dest="structure",
                  help="crystal from file, cif or poscar, REQUIRED",
                  metavar="crystal")

parser.add_option("-r", "--rcut", dest="rcut", default=4.0, type=float,
                  help="cutoff for neighbor calcs, default: 4.0"
                  )

parser.add_option("-l", "--lmax", dest="lmax", default=1, type=int,
                  help="lmax, default: 1"
                  )

parser.add_option("-s", dest="stress", default=True, 
                  action='store_true',help='derivative flag')

parser.add_option("-f", dest="der", default=True,
                  action='store_false',help='derivative flag')

(options, args) = parser.parse_args()

if options.structure is None:
    from ase.build import bulk
    test = bulk('Si', 'diamond', a=5.459)
    cell = test.get_cell()
    cell[0,1] += 0.5
    test.set_cell(cell)
else:
    test = read(options.structure, format='vasp')
print(test)
lmax = options.lmax
rcut = options.rcut
der = options.der
stress = options.stress

#import time
f = SO4_Bispectrum(lmax, rcut, derivative=False, stress=False, normalize_U=False, cutoff_function='tanh')
x = f.calculate(test)
#start2 = time.time()
#for key, item in x.items():
#    print(key, item)
#print('time elapsed: {}'.format(start2 - start1))
#print(x['rdxdr'].shape)
#print(x['rdxdr'])
#print(np.einsum('ijklm->klm', x['rdxdr']))
print(x['x'])
