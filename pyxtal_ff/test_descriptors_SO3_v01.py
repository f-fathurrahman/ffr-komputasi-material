from ase.io import read
import time
from optparse import OptionParser
from my_pyxtal_ff.descriptors_SO3 import SO3

# ---------------------- Options ------------------------
parser = OptionParser()
parser.add_option("-c", "--crystal", dest="structure",
                    help="crystal from file, cif or poscar, REQUIRED",
                    metavar="crystal")

parser.add_option("-r", "--rcut", dest="rcut", default=3.0, type=float,
                    help="cutoff for neighbor calcs, default: 3.0"
                    )

parser.add_option("-l", "--lmax", dest="lmax", default=2, type=int,
                    help="lmax, default: 1"
                    )

parser.add_option("-n", "--nmax", dest="nmax", default=1, type=int,
                    help="nmax, default: 1"
                    )

parser.add_option("-a", "--alpha", dest="alpha", default=2.0, type=float,
                    help="cutoff for neighbor calcs, default: 2.0"
                    )

parser.add_option("-s", dest="stress", default=True,
                    action='store_true',help='derivative flag')

parser.add_option("-f", dest="der", default=True,
                    action='store_false',help='derivative flag')

(options, args) = parser.parse_args()

if options.structure is None:
    from ase.build import bulk
    test = bulk('Si', 'diamond', a=5.459, cubic=True)
    cell = test.get_cell()
    cell[0,1] += 0.5
    test.set_cell(cell)
else:
    test = read(options.structure, format='vasp')

lmax = options.lmax
nmax = options.nmax
rcut = options.rcut
alpha = options.alpha
der = options.der
stress = options.stress

start1 = time.time()
f = SO3(nmax=nmax, lmax=lmax, rcut=rcut, alpha=alpha, derivative=True, stress=False, cutoff_function='cosine')
x = f.calculate(test)
start2 = time.time()
print('x', x['x'])
print('dxdr', x['dxdr'])
print('calculation time {}'.format(start2-start1))
