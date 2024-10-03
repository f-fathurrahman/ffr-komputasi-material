import ase.io
from ase.units import Ry, Bohr
import pickle

fhandle = ase.io.iread("LOG_md01", format="espresso-out")
all_atoms = []
for atoms in fhandle:
    atoms.wrap()
    all_atoms.append(atoms)
    atoms.write("ALL_ATOMS.xyz", append=True)

with open("ALL_ATOMS.pkl", "wb") as f:
    pickle.dump(all_atoms, f)

"""
with open("ALL_ATOMS.pkl", "rb") as f:
    all_atoms = pickle.load(f)
"""
