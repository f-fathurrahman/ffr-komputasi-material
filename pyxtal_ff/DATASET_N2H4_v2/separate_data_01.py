import ase.io

structure_file = "TEMP_ATOMS_TRAIN.xyz"
idata = 1
for atoms in ase.io.iread(structure_file):
    filename = f"STRUC_{idata:04d}.xyz"
    atoms.center(about=[4.0, 4.0, 4.0])
    atoms.write(filename)
    idata += 1
