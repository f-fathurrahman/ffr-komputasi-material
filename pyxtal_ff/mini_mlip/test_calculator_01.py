from optparse import OptionParser
from ase.build import bulk

parser = OptionParser()
parser.add_option("-f", "--file", dest="file",
                    help="pretrained file from pyxtal_ff, REQUIRED",
                    metavar="file")

(options, args) = parser.parse_args()
print(options.file)
calc = MiniMLIPCalculator(mliap=options.file, logo=False)
si = bulk('Si', 'diamond', a=5.459, cubic=True)
si.set_calculator(calc)
print(si.get_potential_energy())
print(si.get_forces())
print(si.get_stress())

box = mushybox(si)
dyn = BFGS(box)
dyn.run(fmax=0.01)
print('equlirum cell para: ', si.get_cell()[0][0])
print('equlirum energy: ', si.get_potential_energy())

Cijs, names, C = elastic_tensor(si, calc)
for name, Cij in zip(names, Cijs):
    print("{:s}: {:8.2f}(GPa)".format(name, Cij))

print(elastic_properties(C))