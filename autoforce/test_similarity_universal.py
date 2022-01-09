import torch
torch.set_default_tensor_type(torch.DoubleTensor)

import numpy as np

from descriptor_atoms import namethem
from descriptor_atoms import TorchAtoms, AtomsData, namethem
from descriptor_cutoff import PolyCut
from similarity_universal import UniversalSoapKernel, DiracDeltaChemical


soap = UniversalSoapKernel(2, 2, 4, 3.)
soap = eval(soap.state)
namethem([soap])

#
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
if 1:
    soap([a], [b]).backward()
    test_left = a.xyz.grad.allclose(soap.leftgrad(a, b).view(-1, 3))
    max_left = (a.xyz.grad - soap.leftgrad(a, b).view(-1, 3)).max()
    print("leftgrad: {}  \t max diff: {}".format(test_left, max_left))
    test_right = b.xyz.grad.allclose(soap.rightgrad(a, b).view(-1, 3))
    max_right = (b.xyz.grad - soap.rightgrad(a, b).view(-1, 3)).max()
    print("rightgrad: {} \t max diff: {}".format(test_right, max_right))

#
if 1:
    data = AtomsData([a, b])
    inducing = data.sample_locals(5)
    inducing.stage(descriptors=[soap])
    soap(data, data)
    soap(data, inducing)
    soap(inducing, inducing)

