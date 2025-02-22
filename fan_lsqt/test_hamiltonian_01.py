import numpy as np
from mod_lsqt_01 import *

Nx = 50000
Ny = 2
W = 1.0

N = Nx * Ny; # total number of sites

row_H = np.zeros(N * 4 - Nx * 2) # the row indices for H
col_H = np.zeros(N * 4 - Nx * 2) # the column indices for H
Hij = -np.ones(N * 4 - Nx * 2, complex) # nonzero Hamiltonian elements

row_V = np.zeros(N * 2) # row indices for V
col_V = np.zeros(N * 2) # column indices for V
Vij = np.zeros(N * 2, complex) # nonzero velocity matrix elements

row_U = np.arange(0, N) # row and column indices for U
Uij = 1.234*np.ones(N) # np.random.uniform(-W * 0.5, W * 0.5, N) # on-site potentials

count_H = 0 # number of nonzero H elements
count_V = 0 # number of nonzero V elements

for nx in range(Nx):
    for ny in range(Ny):
        # (0) # get the index of the center site
        index_center = find_index(nx, ny, Ny)
        # (1) consider the left neighbor (periodic boundary)
        index_left = find_index((nx - 1) % Nx, ny, Ny)
        row_H[count_H] = index_center
        col_H[count_H] = index_left
        count_H += 1
        row_V[count_V] = index_center
        col_V[count_V] = index_left
        Vij[count_V] = 1j
        count_V += 1
        # (2) consider the right neighbor (periodic boundary)
        index_right = find_index((nx + 1) % Nx, ny, Ny)
        row_H[count_H] = index_center
        col_H[count_H] = index_right
        count_H += 1
        row_V[count_V] = index_center
        col_V[count_V] = index_right
        Vij[count_V] = -1j
        count_V += 1
        # (3) consider the upper neighbor (open boundary)
        #print("\n -- nx=%3d ny=%3d index_center=%d" % (nx, ny, index_center))
        #print("Periodic BC: index_left=%3d index_right=%3d" % (index_left, index_right))
        if ny < Ny - 1:
            index_up = find_index(nx, (ny + 1), Ny)
            row_H[count_H] = index_center
            col_H[count_H] = index_up
            count_H += 1
            #print("Open BC: index_up=%3d" % (index_up))
        #else:
        #    print("No index up")
        # (4) consider the down neighbor (open boundary)
        if ny > 0:
            index_down = find_index(nx, (ny - 1), Ny)
            row_H[count_H] = index_center
            col_H[count_H] = index_down
            count_H += 1
            #print("Open BC: index_down=%3d" % (index_down))
        #else:
        #    print("No index down")


H = sparse.csr_matrix((Hij, (row_H, col_H)), shape = (N, N))
U = sparse.csr_matrix((Uij, (row_U, row_U)), shape = (N, N))
H = H + U
V = sparse.csr_matrix((Vij, (row_V, col_V)), shape = (N, N))

#Udense = U.todense()
#Vdense = V.todense()
#Hdense = H.todense()
