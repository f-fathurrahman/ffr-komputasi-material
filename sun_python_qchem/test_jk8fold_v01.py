import numpy as np

from ch12_analytical_integrals_v5.basis import Molecule, num_functions
from ch12_analytical_integrals_v5.eri_OS_rys import get_eri_tensor

from ch13 import j8fold, k8fold, j8fold_v2, k8fold_v2, jk8fold_v2, j_engine


def test_jk8():
    h2o = '''
    O        0.000     0.000     0.000
    H        0.757     0.587     0.000
    H       -0.757     0.587     0.000
    '''
    gtos = Molecule.from_xyz(h2o).assign_basis({'O': '6-31G*', 'H': '6-31G'})

    eri_tensor = get_eri_tensor(gtos)
    nao = eri_tensor.shape[0]
    dm = np.random.rand(nao, nao)
    jref = np.einsum('ijkl,ji->kl', eri_tensor, dm)
    kref = np.einsum('ijkl,jk->il', eri_tensor, dm)

    jmat = j8fold.build_j_8fold(gtos, dm)
    kmat = k8fold.build_k_8fold(gtos, dm)
    assert abs(jref - jmat).max() < 1e-10
    assert abs(kref - kmat).max() < 1e-10
    #
    jmat = j8fold_v2.build_j_8fold(gtos, dm)
    kmat = k8fold_v2.build_k_8fold(gtos, dm)
    assert abs(jref - jmat).max() < 1e-10
    assert abs(kref - kmat).max() < 1e-10
    #
    jmat, kmat = jk8fold_v2.build_jk_8fold(gtos, dm)
    assert abs(jref - jmat).max() < 1e-10
    assert abs(kref - kmat).max() < 1e-10
    #
    print("Pass here in test_jk8")



def test_j_engine():
    h2o = '''
    O        0.000     0.000     0.000
    H        0.757     0.587     0.000
    H       -0.757     0.587     0.000
    '''
    gtos = Molecule.from_xyz(h2o).assign_basis({'O': '6-31G*', 'H': '6-31G'})
    nao = num_functions(gtos)
    dm = np.random.rand(nao, nao)
    jref = j8fold.build_j_8fold(gtos, dm)
    jmat = j_engine.build_j(gtos, dm)
    assert abs(jref - jmat).max() < 1e-10
    print("Pass here in test_j_engine")

test_jk8()
test_j_engine()
