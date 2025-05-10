from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from m3gnet.models import M3GNetCalculator, M3GNet, Potential

# Init a Mo structure with stretched lattice (DFT lattice constant ~ 3.168)
atoms = Structure(
    Lattice.cubic(3.3),
    ["Mo", "Mo"],
    [
        [0., 0., 0.],
        [0.5, 0.5, 0.4]
    ]
)

potential = Potential(M3GNet.load("MP-2021.2.8-EFS"))

# Convert Structure to ase.Atoms ?
atoms = AseAtomsAdaptor().get_atoms(atoms)

stress_weight = 1 / 160.21766208 # what's this?
state_attr = None
element_refs = None
atoms.set_calculator(
    M3GNetCalculator(
        potential=potential,
        stress_weight=stress_weight,
        state_attr=state_attr,
        element_refs=element_refs
    )
)

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
