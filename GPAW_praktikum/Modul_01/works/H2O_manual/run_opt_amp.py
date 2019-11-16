import ase.io
from ase.constraints import FixAtoms
from amp import Amp

atoms = ase.io.read("H2O.xyz")

atoms.center(vacuum=5.0)

atoms.set_pbc([True,True,True])
atoms.write("CENTERED.xsf")

atoms.set_pbc([False,False,False])
atoms.set_constraint(FixAtoms(mask=[0]))

calc = Amp.load("amp.amp")

atoms.set_calculator(calc)

from ase.optimize import BFGS
relax = BFGS( atoms, logfile="LOG_relax_bfgs_AMP", trajectory="geoopt_amp.traj" )
relax.run()
