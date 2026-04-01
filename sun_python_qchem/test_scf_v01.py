from ch12_analytical_integrals_v5.basis import Molecule
from ch13.simple_scf import RHF, scf_iter

xyz = '''
N 0.  0 0
N 1.5 0 0'''
mol = Molecule.from_xyz(xyz)
gtos = mol.assign_basis({'N': '6-31g'})
model = RHF(mol, gtos)
wfn = scf_iter(model)
print('RHF energy', model.total_energy(wfn))

model = RHF.restore(model.chkfile, model.diis.filename)
scf_iter(model, model.wfn)
print('RHF energy', model.total_energy(wfn))

print("Running PySCF")
import pyscf
mol = pyscf.M(atom=list(zip(mol.elements, mol.coordinates)),
                    basis='631g', cart=True, unit='B')
mol.RHF().set(init_guess='hcore').run()
