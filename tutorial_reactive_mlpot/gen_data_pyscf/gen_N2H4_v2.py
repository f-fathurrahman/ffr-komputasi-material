from ase import units
from ase.io import read
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.calculator import Calculator, all_properties
from pyscf import gto, dft, grad


class PySCFCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, restart=None,
                 label='PySCF', atoms=None, scratch=None, **kwargs):
        Calculator.__init__(self, restart, label, atoms, **kwargs)
        self.mol = None
        self.mf = None

    # FIXME: not yet using PBC?
    def calculate(self, atoms=None, properties=['energy', 'forces'],
                  system_changes=all_properties):
        Calculator.calculate(self, atoms, properties, system_changes)
        symbols = atoms.get_chemical_symbols()
        coords = atoms.get_positions()
        self.mol = gto.M(
            atom=[(symbol, tuple(coord)) for symbol, coord in zip(symbols, coords)],
            basis='6-31g(d)'
        )
        self.mf = dft.RKS(self.mol)
        self.mf.xc = 'b3lyp'  # Use B3LYP functional
        self.mf.kernel()
        self.results['energy'] = self.mf.e_tot
        self.results['forces'] = -grad.RKS(self.mf).kernel()



def save_md_results(atoms=None):
    atoms.write("TEMP_langevin_7mols.xyz", append=True)

# Load molecule from file
#atoms = read("N2H4_v1.xyz")
atoms = read('N2H4_7mols.xyz')

# Attach PySCF calculator
atoms.calc = PySCFCalculator()


# Langevin dynamics
T_K = 300
MaxwellBoltzmannDistribution(atoms, temperature_K=T_K)
dyn = Langevin(
    atoms,
    timestep=1.0*units.fs,
    temperature_K=T_K,
    friction=0.001/units.fs,
    logfile="LOG_langevin_7mols"
)
dyn.attach(save_md_results, interval=1, atoms=atoms)
dyn.run(500)
