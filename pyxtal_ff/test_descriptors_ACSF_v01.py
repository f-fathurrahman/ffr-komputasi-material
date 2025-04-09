from ase.build import bulk

from my_pyxtal_ff.descriptors_ACSF import ACSF

import numpy as np
np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

# Set up symmetry parameters
Rc = 5.5
symmetry = {'G2': {'eta': [0.036, 0.071,], 'Rs': [0]},
            'G4': {'Rs': [0], 'lambda': [1], 'zeta': [1,], 'eta': [0.036, 0.071]},
            'G5': {'Rs': [0], 'lambda': [1], 'zeta': [1,], 'eta': [0.036, 0.071]}
           }

for a in [5.0]: #, 5.4, 5.8]:
    si = bulk('Si', 'diamond', a=a, cubic=True)
    cell = si.get_cell()
    cell[0,1] += 0.2
    si.set_cell(cell)
    print(si.get_cell())

    bp = ACSF(symmetry, Rc=Rc, derivative=True, stress=True, cutoff='cosine', atom_weighted=True)
    des = bp.calculate(si, system=[14])
    print(des['x'].shape)
    print("G:", des['x'][0])
    print("Sequence", des['seq'][0])
    print("GPrime", des['dxdr'].shape)
