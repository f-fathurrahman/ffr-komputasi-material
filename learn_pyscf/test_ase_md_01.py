from ase import io, units
from ase.atoms import Atoms
from pyscf_calc import PySCF
from ase.md.verlet import VelocityVerlet
from pyscf import gto, scf

atoms = io.read("N2H4.xyz")
atoms.set_calculator(
    PySCF(atoms=atoms, molcell=gto.M(verbose=0), mf_class=scf.RHF, mf_dict={})
)

dyn = VelocityVerlet(
    atoms, timestep=5.0 * units.fs,
    logfile="TEMP_md.log"
)
dyn.run(100)  # take 1000 steps