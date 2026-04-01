import numpy as np

from ch12_analytical_integrals_v1.basis import CGTO
from ch12_analytical_integrals_v1 import overlap_MD, overlap_OS, overlap_quadrature

prim_gtos = '''
- angular_momentum: 0
  exponents:
  - 0.75
  coefficients:
  - 1.
  coordinates:
  - 0.3
  - 0.0
  - 1.0
- angular_momentum: 1
  exponents:
  - 0.32
  coefficients:
  - 1.
  coordinates:
  - 0.2
  - 1.9
  - 0.0
- angular_momentum: 2
  exponents:
  - 0.51
  coefficients:
  - 1.
  coordinates:
  - 0.5
  - 0.0
  - 0.5
'''

def new_gtos():
    gtos = (CGTO.from_bse('STO-3G', 'C', np.array([0.1, 0.5, 0.8])) +
            CGTO.from_yaml(prim_gtos))
    return gtos

def test_overlap_MD():
    gtos = new_gtos()
    assert abs(abs(overlap_MD.overlap_matrix(gtos)).sum() - 42.91905819232551) < 1e-9

def test_overlap_OS():
    gtos = new_gtos()
    assert abs(abs(overlap_OS.overlap_matrix(gtos)).sum() - 42.91905819232551) < 1e-9

def test_overlap_quadrature():
    gtos = new_gtos()
    assert abs(abs(overlap_quadrature.overlap_matrix(gtos)).sum() - 42.91905819232551) < 1e-9

def test_E_tensor():
    gtos = new_gtos()[3:]
    s = 0
    for bas_i in gtos:
        for bas_j in gtos:
            Et = overlap_MD.get_E_tensor(bas_i.angular_momentum, bas_j.angular_momentum,
                                         bas_i.exponents[0], bas_j.exponents[0],
                                         bas_i.coordinates, bas_j.coordinates)
            s += abs(Et).sum()
    assert abs(s - 53.37624716189072) < 1e-9


from my_utilities_01 import do_timing

do_timing(new_gtos)
do_timing(test_overlap_MD)
do_timing(test_overlap_OS)
do_timing(test_overlap_quadrature)
do_timing(test_E_tensor)
