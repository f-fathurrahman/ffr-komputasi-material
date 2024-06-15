from pyscf import scf as mol_scf
from pyscf.pbc import dft, gto

cell = gto.Cell()
cell.unit = 'Bohr'
cell.atom = '''
Ba  7.643372194672   7.626085433619   7.660733089681
Ti  3.720284319458   3.711870276093   3.728734456591
O   0.165717548281   3.896783227711   3.906001500719
O   3.905615561953   3.896783208814   0.174579400045
O   3.897149929066   0.156894549287   3.906001500719
'''
cell.a = [
[7.622764278, 0.002230425, 0.019531926],
[0.019487064, 7.622738445, 0.019531907],
[0.002237247, 0.002230444, 7.622788202]
]
cell.basis = 'gth-dzvp-molopt-sr'
cell.pseudo = 'gth-pbe'
cell.verbose = 4
cell.build()

kmesh = [2, 2, 2]
kpts = cell.make_kpts(kmesh)

kmf = dft.KRKS(cell, kpts) #.density_fit().apply(mol_scf.addons.remove_linear_dep_)
#kmf.with_df.auxbasis = "weigend"
kmf.xc = 'pbesol'
#kmf = kmf.newton()
kmf.max_cycle = 100
edft = kmf.kernel()

