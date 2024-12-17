import ase.io
from descriptors_SO3 import SO3

atoms = ase.io.read("DATASET_N2H4_v2/N2H4_2mol_1data.xyz")

# Should be invariant with w.r.t translations
pos_shifted = atoms.positions.copy()
pos_shifted[:,2] = atoms.positions[:,2] + 3.0
atoms.set_positions(pos_shifted)

lmax = 4
nmax = 3
rcut = 3.5
alpha = 2.0

desc_calc = SO3(nmax=nmax, lmax=lmax, rcut=rcut, alpha=alpha, derivative=True, stress=False, cutoff_function='cosine')
x = desc_calc.calculate(atoms)

print(desc_calc)
print(x["x"][0])
print(x["x"][1])
print(x["x"].shape)