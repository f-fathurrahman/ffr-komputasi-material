import numpy as np
import ase.io
from ase import Atom, Atoms

atoms = ase.io.read("DATASET_N2H4_v2/N2H4_2mol_1data.xyz")
atsymbs = atoms.get_chemical_symbols()
positions = atoms.get_positions()
cell = atoms.get_cell()
pbc = atoms.pbc
forces = atoms.get_forces()
energy = atoms.get_potential_energy()

idx_sorted = np.argsort(atsymbs)

atoms_copy = Atoms()
for i in idx_sorted:
    atoms_copy.append(Atom(atsymbs[i]))

atoms_copy.set_positions(positions[idx_sorted])
atoms_copy.set_cell(cell)
atoms_copy.set_pbc(atoms.pbc)

from ase.calculators.singlepoint import SinglePointCalculator
stress = None
magmoms = None
calc = SinglePointCalculator(
    atoms_copy,
    energy=energy, forces=forces[idx_sorted]
)
atoms_copy.calc = calc

