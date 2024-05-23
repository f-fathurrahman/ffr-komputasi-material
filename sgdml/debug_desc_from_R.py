import timeit
from functools import partial

import numpy as np
import scipy as sp

def debug_desc_from_R(desc_object, R, lat_and_inv=None):
    
    print("*** ENTER debug_desc_from_R ....")

    # Add singleton dimension if input is (,3N).
    if R.ndim == 1:
        R = R[None, :]

    M = R.shape[0]
    if M == 1:
        return _from_r(R, lat_and_inv)

    R_desc = np.empty([M, desc_object.dim])
    R_d_desc = np.empty([M, desc_object.dim, 3])

    print("R_desc.shape = ", R_desc.shape)
    print("R_d_desc.shape = ", R_d_desc.shape)

    # Generate descriptor and their Jacobians
    start = timeit.default_timer()

    pool = None
    map_func = map
    max_processes = 1

    for i, r_desc_r_d_desc in enumerate(
        map_func(partial(_from_r, lat_and_inv=lat_and_inv), R)
    ):
        R_desc[i, :], R_d_desc[i, :, :] = r_desc_r_d_desc

    stop = timeit.default_timer()
    print("Elapsed time = ", stop - start)

    print("*** EXIT debug_desc_from_R ....")

    return R_desc, R_d_desc



# Generate descriptor and its Jacobian for one molecular geometry
# in Cartesian coordinates.
def _from_r(r, lat_and_inv=None):

    # Add singleton dimension if input is (,3N).
    if r.ndim == 1:
        r = r[None, :]

    pd = _pdist(r, lat_and_inv)

    r_desc = _r_to_desc(r, pd)
    r_d_desc = _r_to_d_desc(r, pd, lat_and_inv)

    return r_desc, r_d_desc


def _pdist(r, lat_and_inv=None):

    r = r.reshape(-1, 3)
    n_atoms = r.shape[0]

    if lat_and_inv is None:
        pdist = sp.spatial.distance.pdist(r, 'euclidean')
    else:
        pdist = sp.spatial.distance.pdist(
            r, lambda u, v: np.linalg.norm(_pbc_diff(u - v, lat_and_inv))
        )

    tril_idxs = np.tril_indices(n_atoms, k=-1)
    return sp.spatial.distance.squareform(pdist, checks=False)[tril_idxs]



def _r_to_desc(r, pdist):
    # Add singleton dimension if input is (,3N).
    if r.ndim == 1:
        r = r[None, :]

    return 1.0 / pdist


# Generate descriptor Jacobian for a set of atom positions in
# Cartesian coordinates.
# 
# This method can apply the minimum-image convention as periodic
# boundary condition for distances between atoms, given the lattice vectors.
def _r_to_d_desc(r, pdist, lat_and_inv=None):
    r = r.reshape(-1, 3)
    pdiff = r[:, None] - r[None, :]  # pairwise differences ri - rj

    n_atoms = r.shape[0]
    i, j = np.tril_indices(n_atoms, k=-1)

    pdiff = pdiff[i, j, :]  # lower triangular

    if lat_and_inv is not None:
        pdiff = _pbc_diff(pdiff, lat_and_inv)

    d_desc_elem = pdiff / (pdist ** 3)[:, None]

    return d_desc_elem


# Clamp differences of vectors to super cell.
def _pbc_diff(diffs, lat_and_inv):

    lat, lat_inv = lat_and_inv

    c = lat_inv.dot(diffs.T)
    diffs -= lat.dot(np.around(c)).T

    return diffs
