import ase.io

ATOMS_LIST = []
for atoms in ase.io.iread("OUTCAR_comp"):
    ATOMS_LIST.append(atoms)


