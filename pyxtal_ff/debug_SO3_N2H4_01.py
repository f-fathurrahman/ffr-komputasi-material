import ase.io
from my_desc_SO3 import MyDescriptorSO3, calculate_SO3_power_spectrum

atoms = ase.io.read("DATASET_N2H4_v2/N2H4_2mol_1data.xyz")

#atoms_list = ase.io.read("DATASET_OTHERS/TiAl_gabung.xyz@:")
#atoms = atoms_list[0]

# Should be invariant with w.r.t translations
#pos_shifted = atoms.positions.copy()
#pos_shifted[:,2] = atoms.positions[:,2] + 3.0
#atoms.set_positions(pos_shifted)

print("Number of atoms = ", len(atoms))

lmax = 4
nmax = 3
rcut = 3.5
alpha = 2.0

desc_calc = MyDescriptorSO3(
    nmax=nmax,
    lmax=lmax,
    rcut=rcut,
    alpha=alpha,
    derivative=False,
    stress=False,
    cutoff_function='cosine'
)


# desc_calc.center_atoms = np.array(center_atoms, dtype=np.float64)
# desc_calc.neighborlist = np.array(neighbors, dtype=np.float64)
# desc_calc.seq = Seq
# desc_calc.atomic_weights = np.array(atomic_weights, dtype=np.int64)
# desc_calc.neighbor_indices = neighbor_indices

#x = desc_calc.calculate(atoms)
x = calculate_SO3_power_spectrum(desc_calc, atoms)
print(desc_calc)

print(x["x"][0])
print(x["x"][1])
print(x["x"].shape)
#print(x["seq"])
