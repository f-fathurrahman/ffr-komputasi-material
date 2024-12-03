import ase.io

structure_file = "DATASET_N2H4_v1/TEMP_ATOMS_TRAIN.xyz"
data = []
for atoms in ase.io.iread(structure_file):
    data.append({
        'structure': atoms,
        'energy': atoms.get_potential_energy(),
        'force': atoms.get_forces(),
        'stress': None,
        'group': 'random'
    })
