import timeit
from functools import partial

import numpy as np
import scipy as sp

# removed references to Desc object


def my_desc_from_R(desc_dim, R, lat_and_inv=None):

    # desc_dim is Natoms*(Natoms-1) // 2

    # Add singleton dimension if input is (,3N).
    if R.ndim == 1:
        R = R[None, :]

    M = R.shape[0]
    # only one example data
    if M == 1:
        R_desc, R_d_desc = _from_r(R, lat_and_inv)
        print("R_desc.shape = ", R_desc.shape)
        print("R_d_desc.shape = ", R_d_desc.shape)
        return R_desc, R_d_desc

    R_desc = np.empty([M, desc_dim])
    R_d_desc = np.empty([M, desc_dim, 3])

    map_func = map

    for i, r_desc_r_d_desc in enumerate(
        map_func(partial(_from_r, lat_and_inv=lat_and_inv), R)
    ):
        R_desc[i, :], R_d_desc[i, :, :] = r_desc_r_d_desc

    print("R_desc.shape = ", R_desc.shape)
    print("R_d_desc.shape = ", R_d_desc.shape)
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

    # XXX: Why need this?
    tril_idxs = np.tril_indices(n_atoms, k=-1)
    return sp.spatial.distance.squareform(pdist, checks=False)[tril_idxs]



def _r_to_desc(r, pdist):
    # Add singleton dimension if input is (,3N).
    # XXX: Why?
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
    # pdiff is an array of vector differences

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



def d_desc_from_comp(n_atoms, R_d_desc, out=None):

    desc_dim = (n_atoms * (n_atoms - 1)) // 2
    tril_indices = np.tril_indices(n_atoms, k=-1)
    dim_i = 3 * n_atoms

    if R_d_desc.ndim == 2:
        R_d_desc = R_d_desc[None, ...]

    n = R_d_desc.shape[0] # number of data
    i, j = tril_indices

    if out is None:
        out = np.zeros((n, desc_dim, n_atoms, 3))
    else:
        out = out.reshape(n, desc_dim, n_atoms, 3)

    dim_range = np.arange(desc_dim)  # [0, 1, ..., dim-1]
    out[:, dim_range, j, :] = R_d_desc
    out[:, dim_range, i, :] = -R_d_desc

    return out.reshape(-1, desc_dim, dim_i)

