from ase.build import bulk

Cu_v1 = bulk("Cu", "fcc", a=3.6)
print(Cu_v1.cell)
Cu_v1.write("Cu_v1.xsf")

Cu_v2 = bulk("Cu", "fcc", a=3.6, orthorhombic=True)
print(Cu_v2.cell)
Cu_v2.write("Cu_v2.xsf")

Cu_v3 = bulk("Cu", "fcc", a=3.6, cubic=True)
print(Cu_v3.cell)
Cu_v3.write("Cu_v3.xsf")