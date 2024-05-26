from my_sgdml.intf.ase_calc import SGDMLCalculator

from ase.io import read
from ase.optimize import QuasiNewton
from ase.vibrations import Vibrations

#model_path = "PRETRAINED/ethanol_1k.npz"
model_path = "PRETRAINED/ethanol_50k.npz" # take much longer
#model_path = "m_ethanol.npz"
calc = SGDMLCalculator(model_path)

mol = read("ethanol.xyz")
mol.set_calculator(calc)

# do a quick geometry relaxation
qn = QuasiNewton(mol)
qn.run(1e-4, 100)

# run the vibration calculations
vib = Vibrations(mol)
vib.run()

vib.summary() # print a summary of the vibrational frequencies
vib.write_jmol() # write file for viewing of the modes with jmol

vib.clean() # remove pickle-files