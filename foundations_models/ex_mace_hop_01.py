from ase.build import add_adsorbate, fcc100, fcc111, fcc110, fcc211, molecule
from ase.constraints import FixAtoms
from ase.optimize import LBFGS
from mace.calculators import MACECalculator
from ase.constraints import FixAtoms, Hookean
from ase.optimize.minimahopping import MinimaHopping

atsymb = "Ni"
slab = fcc110(atsymb, (2, 3, 5), vacuum=8, periodic=True)

# Add adsorbate at the center
adsorbate = molecule("N2H4")
cell = slab.get_cell()
xc, yc, _ = 0.5 * (cell[0] + cell[1])
height = 2.0
add_adsorbate(slab, adsorbate, height, position=(xc - 1.0, yc))

slab.write("GEOM.xyz")
#from ase.visualize import view
#view(slab)
#atoms.extend(adsorbate) # Need this?

# N-N bond in hydrazine is about 1.46 angstrom (from wikipedia)
zmin = slab.get_positions()[:, 2].min()
#
bond_N1_N2 = Hookean(a1=30, a2=31, rt=1.8, k=15.0)
#
#bond_N1_H1 = Hookean(a1=30, a2=32, rt=1.2, k=15.0)
bond_N1_H2 = Hookean(a1=30, a2=33, rt=1.2, k=15.0)
#
bond_N2_H3 = Hookean(a1=31, a2=34, rt=1.2, k=15.0)
bond_N2_H4 = Hookean(a1=31, a2=35, rt=1.2, k=15.0)
#
constraints = [
    FixAtoms(mask=slab.positions[:, 2] < zmin + 0.01),
    bond_N1_N2,
    bond_N1_H2,
    bond_N2_H3, bond_N2_H4,
]
slab.set_constraint(constraints)

# Set the calculator.
model_path = "../pretrained/mace-mh-1.model"
calc = MACECalculator(
    model_paths="../pretrained/mace-mh-1.model", device="cpu", head="oc20_usemppbe"
)

slab.calc = calc

# Instantiate and run the minima hopping algorithm.
hop = MinimaHopping(slab, Ediff0=2.5, T0=4000.0)
#dyn.attach(atoms.wrap, interval=1) # XXX probably need this for long totalsteps
hop(totalsteps=10)

from ase.optimize.minimahopping import MHPlot
mhplot = MHPlot()
mhplot.save_figure('summary.png')
