from ase.build import bulk

atoms = bulk("Ni", cubic=True, a=4.52)
atoms.write("Ni.xyz")
atoms.write("Ni.xsf")
