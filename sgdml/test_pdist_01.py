import numpy as np
from scipy.spatial.distance import pdist

Natoms = 5
R = np.zeros((Natoms,3))
R[1] = [0.0, 0.0, 1.1]
R[2] = [0.0, 0.0,  0.1]
R[3] = [1.0, 0.0, 90.0]
R[4] = [0.0, 0.0, 11.0]

# Total number of elements of dR: Natoms*(Natoms-1)/2
dR = pdist(R)
""""
Returns:
Y ndarray
Returns a condensed distance matrix Y. For each i and j (where i < j < m),
where m is the number of original observations.
The metric dist(u=X[i], v=X[j]) is computed and stored in entry m * i + j - ((i + 2) * (i + 1)) // 2.
"""

dR_manual = np.zeros( int(Natoms*(Natoms-1)/2) )
ip = 0
for i in range(Natoms):
    for j in range(i+1,Natoms):
        print(f"{i} {j}")
        dR_manual[ip] = np.linalg.norm(R[i] - R[j])
        ip += 1
