import ase.io
from gpaw import GPAW, PW, Davidson, FermiDirac
from ase.constraints import FixAtoms

atoms = ase.io.read("H2O.xyz")

atoms.center(vacuum=5.0)

atoms.set_pbc([True,True,True]) # set temporarily to True for visualization
atoms.write("CENTERED.xsf")

atoms.set_pbc([False,False,False])
atoms.set_constraint(FixAtoms(mask=[0]))

# gpaw calculator:
calc = GPAW( mode="lcao", basis="dzp",
             #mode=PW(ecut=340),
             xc="PBE",
             occupations=FermiDirac(0.001),
             txt="-" )

atoms.set_calculator(calc)

#atoms.get_forces()
#calc.write("restart.gpw")

from ase.optimize import BFGS, QuasiNewton

relax = BFGS( atoms, logfile="LOG_relax_bfgs_LCAO", trajectory="geoopt_LCAO.traj" )
#relax = BFGS( atoms, logfile="LOG_relax_bfgs_PW", trajectory="geoopt_PW.traj" )

relax.run()
#relax.run(fmax=0.05)

calc.write("restart_opt.gpw")