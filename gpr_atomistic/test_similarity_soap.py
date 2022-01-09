import torch

from similarity_soap import SoapKernel, NormedSoapKernel

def test_grad():
    from descriptor_atoms import namethem
    from descriptor_cutoff import PolyCut
    from regression_kernel import Positive, DotProd, Normed
    from regression_stationary import RBF
    from descriptor_atoms import TorchAtoms, AtomsData
    import numpy as np
    torch.set_default_tensor_type(torch.DoubleTensor)

    # create kernel
    kern = Positive(1.0) * Normed(DotProd())**4
    #kern = RBF()
    #soap = SoapKernel(kern, 10, (18, 10), 2, 2, PolyCut(3.0))
    soap = NormedSoapKernel(DotProd()**4, 10, (18, 10), 2, 2, PolyCut(3.0))
    namethem([soap])

    # create atomic systems
    # Note that when one of the displacement vectors becomes is exactly along the z-axis
    # because of singularity some inconsistensies exist with autograd.
    # For this reason we add a small random number to positions, until that bug is fixed.
    cell = np.ones(3)*10
    positions = np.array([(-1., 0., 0.), (1., 0., 0.),
                          (0., -1., 0.), (0., 1., 0.),
                          (0., 0., -1.), (0., 0., 1.0),
                          (0., 0., 0.)]) + cell/2 + np.random.uniform(-0.1, 0.1, size=(7, 3))

    b = TorchAtoms(positions=positions, numbers=3*[10]+3*[18]+[10], cell=cell,
                   pbc=True, cutoff=3.0, descriptors=[soap])
    # make natoms different in a, b. P.S. add an isolated atom.
    _pos = np.concatenate([positions, [[0., 0., 0.], [3., 5., 5.]]])
    a = TorchAtoms(positions=_pos, numbers=2*[10, 18, 10]+[18, 10, 18], cell=cell,
                   pbc=True, cutoff=3.0, descriptors=[soap])

    # left/right-grad
    a.update(posgrad=True, forced=True)
    b.update(posgrad=True, forced=True)
    soap([a], [b]).backward()
    test_left = a.xyz.grad.allclose(soap.leftgrad(a, b).view(-1, 3))
    max_left = (a.xyz.grad - soap.leftgrad(a, b).view(-1, 3)).max()
    print("leftgrad: {}  \t max diff: {}".format(test_left, max_left))
    test_right = b.xyz.grad.allclose(soap.rightgrad(a, b).view(-1, 3))
    max_right = (b.xyz.grad - soap.rightgrad(a, b).view(-1, 3)).max()
    print("rightgrad: {} \t max diff: {}".format(test_right, max_right))

    # gradgrad-left
    a.update(posgrad=True, forced=True)
    b.update(posgrad=True, forced=True)
    (soap.leftgrad(a, b).view(-1, 3)*a.xyz).sum().backward()
    v1 = a.xyz.grad.data
    a.update(posgrad=True, forced=True)
    b.update(posgrad=True, forced=True)
    (soap.gradgrad(a, b)*a.xyz.view(-1)[:, None]).sum().backward()
    v2 = a.xyz.grad.data
    print('gradgrad-left: {}'.format(v1.allclose(v2)))

    # gradgrad-right
    a.update(posgrad=True, forced=True)
    b.update(posgrad=True, forced=True)
    (soap.rightgrad(a, b).view(-1, 3)*b.xyz).sum().backward()
    v1 = b.xyz.grad.data
    a.update(posgrad=True, forced=True)
    b.update(posgrad=True, forced=True)
    (soap.gradgrad(a, b)*b.xyz.view(-1)[None]).sum().backward()
    v2 = b.xyz.grad.data
    print('gradgradright: {}'.format(v1.allclose(v2)))

    # gradgraddiag
    test_diag = soap.gradgrad(a, a).diag().allclose(soap.gradgraddiag(a))
    print('gradgraddiag: {}'.format(test_diag))


def example():
    from regression_kernel import Positive, DotProd, Mul, Add, Pow
    from descriptor_cutoff import PolyCut
    kern = (Positive(1.0, requires_grad=True) *
            (DotProd() + Positive(0.01, requires_grad=True))**0.1)
    soap = SoapKernel(kern, 10, (18, 10), 2, 2, PolyCut(3.0))
    assert eval(soap.state).state == soap.state


example()
test_grad()
