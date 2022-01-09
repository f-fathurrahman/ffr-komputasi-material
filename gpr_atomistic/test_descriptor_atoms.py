import numpy as np
from similarity_pair import DistanceKernel
from regression_core import SquaredExp
from descriptor_atoms import TorchAtoms, namethem

kerns = [DistanceKernel(SquaredExp(), 10, 10),
         DistanceKernel(SquaredExp(), 10, 18),
         DistanceKernel(SquaredExp(), 18, 18)]
namethem(kerns)
xyz = np.stack(np.meshgrid([0, 1.5], [0, 1.5], [0, 1.5])
               ).reshape(3, -1).transpose()
numbers = 4*[10] + 4*[18]
atoms = TorchAtoms(positions=xyz, numbers=numbers,
                   cutoff=3.0, descriptors=kerns)

other = atoms.copy()
print(other == atoms)

for loc in atoms:
    print(loc.as_atoms().as_local() == loc.detach())

empty = TorchAtoms(positions=[(0, 0, 0)], cutoff=3.)
empty[0].detach()._r
