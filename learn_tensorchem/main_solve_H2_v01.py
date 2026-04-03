import numpy as np
from qtools.molecules import Molecule
import tucker3d as tuck
from tensorchem_solvers.solver import Solver

name = "H2"
mol = Molecule(name)

n = 2**11; print("n = ", n)
sol = Solver(mol, method='hf', mixing=2, eps=1e-6, maxiter=60, meshsize=n, boxsize=12.0)
sol.solve()