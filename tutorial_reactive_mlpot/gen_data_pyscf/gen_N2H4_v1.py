from ase import units, Atoms
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.calculator import Calculator, all_properties
from ase.optimize import QuasiNewton
import numpy as np
from pyscf import gto, dft, grad


class PySCFCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, restart=None,
                 label='PySCF', atoms=None, scratch=None, **kwargs):
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)
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




# Load molecule from file
atoms = read('N2H4_7mols.xyz')

# Menentukan unit cell dengan panjang sisi 20 Å
unit_cell_size = 20.0
atoms.set_cell(unit_cell_size * np.identity(3))
atoms.set_pbc(True)  # Mengatur Periodic Boundary Conditions (PBC)

# Hitung pusat massa (center of mass)
center_of_mass = atoms.get_center_of_mass()

# Hitung vektor translasi agar pusat massa berada di tengah unit cell
translation_vector = 0.5 * unit_cell_size - center_of_mass

# Pindahkan semua atom ke posisi baru berdasarkan vektor translasi
atoms.translate(translation_vector)

# Pastikan semua atom berada di dalam unit cell setelah translasi
atoms.wrap()

# Attach PySCF calculator
atoms.calc = PySCFCalculator()

atoms.get_total_energy()


