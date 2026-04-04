from qtools.molecules import Molecule
from tensorchem_solvers.solver import Solver

name = "H2O"
mol = Molecule(name)

n = 2**11; print("n = ", n)
#sol = Solver(mol, method='hf', mixing=2, eps=1e-6, maxiter=60, meshsize=n, boxsize=12.0)
sol = Solver(mol, method='lda', mixing=2, eps=1e-6, maxiter=60, meshsize=n, boxsize=12.0)
sol.solve()