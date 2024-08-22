# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 2
#     language: python
#     name: python2
# ---

# %% [markdown]
# # Basis Sets 

# %% run_control={"frozen": false, "read_only": false}
import psi4
import numpy as np

# %% run_control={"frozen": false, "read_only": false}
# ==> Set Basic Psi4 Options <==
# Memory specification
psi4.set_memory(int(5e8))
numpy_memory = 2

# Set output file
psi4.core.set_output_file('output.dat', False)

# %% [markdown] run_control={"frozen": false, "read_only": false}
# Thus far you've used a uniform single basis set from [Psi4's basis set library](https://github.com/psi4/psi4/tree/master/psi4/share/psi4/basis) where you've specified the orbital basis (*e.g.*, `psi4.set_options({'basis': '6-31G*'})` and allowed Psi4 to select any necessary auxiliary basis sets to fulfill your requested algorithm. In this tutorial, we'll learn how to construct custom orbital and auxiliary basis set and get information from the psi4.core.BasisSet object

# %% [markdown] run_control={"frozen": false, "read_only": false}
# One distinction that's important to make early on is that a "BasisSet object" is always tailored to a Molecule --- there are no shells assigned to carbon in a cc-pVDZ BasisSet associated with water. In contrast, a "basis set definition" is a set of rules for applying shells to many Molecules (not _any_ Molecule-s because the definition mightn't include uranium, for instance) on the basis of elemental composition and atom labeling. There's nothing stopping you from assigning carbon-parametrized shells to carbon _and_ oxygen in a basis set definition. When the basis set definition is applied to water, the resulting BasisSet object will have carbon-parametrized shells assigned to oxygen but no shells assigned to carbon. Keep this distinction in mind since a basis set like `cc-pVDZ` is commonly used in both roles interchangeably.

# %% run_control={"frozen": false, "read_only": false}
from pkg_resources import parse_version
if parse_version(psi4.__version__) >= parse_version('1.3a1'):
    refnuc =  204.01995818060678
    refscf = -228.95763005900784
else:
    refnuc =  204.01995737868003
    refscf = -228.95763005849557

bzb = psi4.geometry("""
    X
    X   1  RXX
    X   2  RXX  1  90.0
    C   3  RCC  2  90.0  1   0.0
    C   3  RCC  2  90.0  1  60.0
    C1  3  RCC  2  90.0  1 120.0
    C   3  RCC  2  90.0  1 180.0
    C1  3  RCC  2  90.0  1 240.0
    C   3  RCC  2  90.0  1 300.0
    H1  3  RCH  2  90.0  1   0.0
    H   3  RCH  2  90.0  1  60.0
    H   3  RCH  2  90.0  1 120.0
    H1  3  RCH  2  90.0  1 180.0
    H   3  RCH  2  90.0  1 240.0
    H   3  RCH  2  90.0  1 300.0
    RCC  = 1.3915
    RCH  = 2.4715
    RXX  = 1.00
""")
psi4.set_options({'mp_type': 'conv'})
psi4.core.IO.set_default_namespace("bzb")

def basisspec_psi4_yo__anonymous775(mol, role):
    basstrings = {}
    mol.set_basis_all_atoms("DZ", role=role)
    mol.set_basis_by_symbol("C", "my3-21G", role=role)
    mol.set_basis_by_label("H1", "sto-3g", role=role)
    mol.set_basis_by_label("C1", "sto-3g", role=role)
    basstrings['my3-21g'] = """
cartesian
****
H     0
S   2   1.00
5.4471780              0.1562850
0.8245470              0.9046910
S   1   1.00
0.1831920              1.0000000
****
C     0
S   3   1.00
172.2560000              0.0617669
25.9109000              0.3587940
5.5333500              0.7007130
SP   2   1.00
3.6649800             -0.3958970              0.2364600
0.7705450              1.2158400              0.8606190
SP   1   1.00
0.1958570              1.0000000              1.0000000
****
"""
    basstrings['dz'] = """
spherical
****
H     0
S   3   1.00
19.2406000              0.0328280
2.8992000              0.2312080
0.6534000              0.8172380
S   1   1.00
0.1776000              1.0000000
****
"""
    return basstrings

psi4.qcdb.libmintsbasisset.basishorde['ANONYMOUS775'] = basisspec_psi4_yo__anonymous775

psi4.set_options({'basis': 'anonymous775',
                  'scf_type': 'pk',
                  'e_convergence': 11,
                  'd_convergence': 11})

eb, wb = psi4.energy('scf', return_wfn=True)
psi4.compare_strings("c2v", bzb.schoenflies_symbol(), "Point group")                       
psi4.compare_values(refnuc, bzb.nuclear_repulsion_energy(), 10, "Nuclear repulsion energy")
psi4.compare_values(refscf, eb, 7, "SCF Energy")                                  

# %% run_control={"frozen": false, "read_only": false}
psi4.core.set_output_file('output.dat', True)

bsb = wb.get_basisset('ORBITAL')
bsb.print_detail_out()
bsb.print_out()

# %% run_control={"frozen": false, "read_only": false}
#           cc-pvdz                 aug-cc-pvdz
# BASIS     H  5/ 5   C  14/15      H +4/ 4   C  +9/10
# RIFIT     H 14/15   C  56/66      H +9/10   C +16/20
# JKFIT     H 23/25   C  70/81      H +9/10   C +16/20

mymol = psi4.qcdb.Molecule("""
C    0.0  0.0 0.0
O    1.4  0.0 0.0
H_r -0.5 -0.7 0.0
H_l -0.5  0.7 0.0
""")

print('[1]    <<<  uniform cc-pVDZ  >>>')
wert = psi4.qcdb.BasisSet.pyconstruct(mymol, 'BASIS', 'cc-pvdz')
psi4.compare_integers(38, wert.nbf(), 'nbf()')
psi4.compare_integers(40, wert.nao(), 'nao()')
psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')
mymol.print_out()


print('[2]        <<<  RIFIT (default)  >>>')
wert = psi4.qcdb.BasisSet.pyconstruct(mymol, 'DF_BASIS_MP2', '', 'RIFIT', 'cc-pvdz')
psi4.compare_integers(140, wert.nbf(), 'nbf()')
psi4.compare_integers(162, wert.nao(), 'nao()')
psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')
mymol.print_out()


print('[3]    <<<  cc-pVDZ w/ aug-cc-pVDZ on C  >>>')
def basisspec_psi4_yo__anonymous775(mol, role):
    basstrings = {}
    mol.set_basis_all_atoms("DZ", role=role)
    mol.set_basis_by_symbol("C", "my3-21G", role=role)

def basis__dz_PLUS(mol, role):
    mol.set_basis_all_atoms("cc-pvdz", role=role)
    mol.set_basis_by_symbol("C", "aug-cc-pvdz")
    return {}

wert = psi4.qcdb.BasisSet.pyconstruct(mymol, 'BASIS', basis__dz_PLUS)
psi4.compare_integers(47, wert.nbf(), 'nbf()')
psi4.compare_integers(50, wert.nao(), 'nao()')
psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')
mymol.print_out()

