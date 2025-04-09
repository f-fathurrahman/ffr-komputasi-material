import numpy as np
from ase.build import bulk
from my_pyxtal_ff.descriptors_EAD import EAD

np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

Rc = 10
parameters1 = {'L': 2, 'eta': [0.036, 0.071], 'Rs': [0]}

# Test for stress
for a in [5.0]: #, 5.4, 5.8]:
    si = bulk('Si', 'diamond', a=a, cubic=True)
    cell = si.get_cell()
    cell[0,1] += 0.1
    si.set_cell(cell)

    bp = EAD(parameters1, Rc=Rc, derivative=True, stress=True, cutoff='cosine')
    des = bp.calculate(si)
    
    print("G:", des['x'])
    print("GPrime")
    print(des['dxdr'])
    #print(des['rdxdr'][0:8, -1, :, :])
    #pprint(np.einsum('ijklm->klm', des['rdxdr']))
