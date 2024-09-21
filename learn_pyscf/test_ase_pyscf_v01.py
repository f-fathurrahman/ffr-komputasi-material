from ase import io
from ase.atoms import Atoms
from pyscf_simple import PySCF_simple

mol = io.read("N2H4.xyz")
mol.set_calculator(PySCF_simple(atoms=mol, method='DFT', basis='6-31g*'))

print(mol.get_potential_energy())
print(mol.get_forces())

