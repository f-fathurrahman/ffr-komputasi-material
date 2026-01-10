from ase.build import add_adsorbate, fcc100, fcc111, fcc211, molecule
from ase.constraints import FixAtoms
from ase.optimize import LBFGS
from fairchem.core import FAIRChemCalculator
from fairchem.core.units.mlip_unit import load_predict_unit

# Prepare calculator
model_path = "../uma-s-1p1.pt"
predictor = load_predict_unit(path=model_path, device="cpu")
calc = FAIRChemCalculator(predictor, task_name="oc20")


def run_optim(atsymb, orientation, sx, sy, itry):
    # periodic=True (for all directions) is needed for FAIRChemCalculator
    if orientation == "111":
        slab = fcc111(atsymb, (3, 3, 3), vacuum=8, periodic=True)
    elif orientation == "100":
        slab = fcc100(atsymb, (3, 3, 3), vacuum=8, periodic=True)
    elif orientation == "211":
        slab = fcc211(
            atsymb, (3, 3, 3), vacuum=8
        )  # XXX: kwarg periodic is not supported?
        slab.set_pbc([True, True, True])
    else:
        raise RuntimeError(f"Not supported: {atsymb} for orientation {orientation}")

    adsorbate = molecule("N2H4")

    # Add adsorbate at the center
    cell = slab.get_cell()
    xc, yc, _ = 0.5 * (cell[0] + cell[1])
    height = 2.0
    add_adsorbate(slab, adsorbate, height, position=(xc + sx, yc + sy))
    # Fix lowest layer of the slab
    zmin = slab.get_positions()[:, 2].min()
    constraint = FixAtoms(mask=slab.positions[:, 2] < zmin + 0.01)
    slab.set_constraint(constraint)

    slab.calc = calc
    # Set up LBFGS dynamics object
    ftraj = f"{atsymb}{orientation}_N2H4_{itry}.traj"
    opt = LBFGS(slab, trajectory=ftraj)
    opt.run(0.05, 1000)


atsymb = "Ni"  # Ni, Cu, Pt, Rh, Ir, Pd
orientation = "100"
itry = 0
for sx in [-1, 0, 1]:
    for sy in [-1, 0, 1]:
        itry += 1
        run_optim(atsymb, orientation, sx, sy, itry)
