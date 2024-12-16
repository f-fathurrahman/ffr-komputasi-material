import ase.io

structure_file = "TiAl_gabung.xyz"
idata = 1
for atoms in ase.io.iread(structure_file):
    filename = f"STRUC_{idata:04d}.xyz"
    atoms.write(filename)
    idata += 1
