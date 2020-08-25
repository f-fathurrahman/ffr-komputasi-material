import copy
from ase import Atoms
from MullerBrown import MullerBrown

from ase.build import fcc100, add_adsorbate
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms

from do_optim_BFGS import do_optim_BFGS

def prepare_MullerBrown():
    calc = MullerBrown()
    # "Molecular" structures:
    # Using one atom only, position at xy plane only (with z=0)
    atoms = Atoms("C", positions=[(-0.55, 1.30, 0.0)])
    atoms.set_calculator(copy.deepcopy(calc))
    #
    return atoms

def prepare_Au_on_Al_slab():
    atoms = fcc100("Al", size=(2, 2, 3))
    add_adsorbate(atoms, "Au", 1.7, "hollow")
    atoms.center(axis=2, vacuum=4.0)
    atoms.set_calculator(EMT())
    # Fix second and third layers:
    mask = [atom.tag > 1 for atom in atoms]
    atoms.set_constraint(FixAtoms(mask=mask))
    #
    return atoms


#atoms = prepare_MullerBrown()
atoms = prepare_Au_on_Al_slab()

# Test calculating potential energy and forces
print("Potential energy = ", atoms.get_potential_energy())
print("Forces           = ")
print(atoms.get_forces())

print("Optimization\n")
print("Initial positions: ")
print(atoms.get_positions())

do_optim_BFGS(atoms)

print("Optimized positions: ")
print(atoms.get_positions())
