import ase.io

atoms = ase.io.read("Cu.cif")
Natoms = len(atoms)
atpos = atoms.positions
cell = atoms.cell

print("ATOMIC_POSITIONS bohr")
for ia in range(Natoms):
    print("%s %18.10f %18.10f %18.10f" % (atoms[ia].symbol, atpos[ia,0], atpos[ia,1], atpos[ia,2]))

print("")
print("CELL_PARAMETERS bohr")
for i in range(3):
    print("%18.10f %18.10f %18.10f" % (cell[i,0], cell[i,1], cell[i,2]))

