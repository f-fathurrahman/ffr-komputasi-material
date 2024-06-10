import numpy as np
import scipy as sp
from functools import partial


def _pbc_diff(diffs, lat_and_inv):
    lat, lat_inv = lat_and_inv
    c = lat_inv.dot(diffs.T)
    diffs -= lat.dot(np.around(c)).T
    return diffs


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


def _squareform(vec_or_mat):

    # vector to matrix representation
    if vec_or_mat.ndim == 1:

        n_tril = vec_or_mat.size
        n = int((1 + np.sqrt(8 * n_tril + 1)) / 2)

        i, j = np.tril_indices(n, k=-1)

        mat = np.zeros((n, n))
        mat[i, j] = vec_or_mat
        mat[j, i] = vec_or_mat

        return mat

    else:  # matrix to vector

        assert vec_or_mat.shape[0] == vec_or_mat.shape[1]  # matrix is square

        n = vec_or_mat.shape[0]
        i, j = np.tril_indices(n, k=-1)

        return vec_or_mat[i, j]


def _r_to_desc(r, pdist):
    # Add singleton dimension if input is (,3N).
    if r.ndim == 1:
        r = r[None, :]

    return 1.0 / pdist


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


def _from_r(r, lat_and_inv=None):

    # Add singleton dimension if input is (,3N).
    if r.ndim == 1:
        r = r[None, :]

    pd = _pdist(r, lat_and_inv)

    r_desc = _r_to_desc(r, pd)
    r_d_desc = _r_to_d_desc(r, pd, lat_and_inv)

    return r_desc, r_d_desc


class Desc(object):

    def __init__(self, n_atoms, max_processes=None):
        self.n_atoms = n_atoms
        self.dim_i = 3 * n_atoms

        # Size of the resulting descriptor vector.
        self.dim = (n_atoms * (n_atoms - 1)) // 2

        self.tril_indices = np.tril_indices(n_atoms, k=-1)

        # Precompute indices for nonzero entries in desriptor derivatives.
        self.d_desc_mask = np.zeros((n_atoms, n_atoms - 1), dtype=np.int32)
        for a in range(n_atoms):  # for each partial derivative
            rows, cols = self.tril_indices
            self.d_desc_mask[a, :] = np.concatenate(
                [np.where(rows == a)[0], np.where(cols == a)[0]]
            )

        self.dim_range = np.arange(self.dim)  # [0, 1, ..., dim-1]

        # Precompute indices for nonzero entries in desriptor derivatives.

        self.M = np.arange(1, n_atoms)  # indexes matrix row-wise, skipping diagonal
        for a in range(1, n_atoms):
            self.M = np.concatenate((self.M, np.delete(np.arange(n_atoms), a)))

        self.A = np.repeat(
            np.arange(n_atoms), n_atoms - 1
        )  # [0, 0, ..., 1, 1, ..., 2, 2, ...]

        self.max_processes = max_processes

    def from_R(self, R, lat_and_inv=None, max_processes=None, callback=None):

        # Add singleton dimension if input is (,3N).
        if R.ndim == 1:
            R = R[None, :]

        M = R.shape[0]
        if M == 1:
            return _from_r(R, lat_and_inv)

        R_desc = np.empty([M, self.dim])
        R_d_desc = np.empty([M, self.dim, 3])

        map_func = map
        for i, r_desc_r_d_desc in enumerate(
            map_func(partial(_from_r, lat_and_inv=lat_and_inv), R)
        ):
            R_desc[i, :], R_d_desc[i, :, :] = r_desc_r_d_desc

        return R_desc, R_d_desc


    # Multiplies descriptor(s) jacobian with 3N-vector(s) from the right side
    def d_desc_dot_vec(self, R_d_desc, vecs, overwrite_vecs=False):

        if R_d_desc.ndim == 2:
            print("Pass here 154 in my_desc")
            R_d_desc = R_d_desc[None, ...]

        if vecs.ndim == 1:
            print("Pass here 156 in my_desc")
            vecs = vecs[None, ...]

        i, j = self.tril_indices

        vecs = vecs.reshape(vecs.shape[0], -1, 3)

        einsum = np.einsum

        return einsum('...ij,...ij->...i', R_d_desc, vecs[:, j, :] - vecs[:, i, :])

    # Multiplies descriptor(s) jacobian with N(N-1)/2-vector(s) from the left side
    def vec_dot_d_desc(self, R_d_desc, vecs, out=None):

        if R_d_desc.ndim == 2:
            R_d_desc = R_d_desc[None, ...]

        if vecs.ndim == 1:
            vecs = vecs[None, ...]

        assert (
            R_d_desc.shape[0] == 1
            or vecs.shape[0] == 1
            or R_d_desc.shape[0] == vecs.shape[0]
        )  # either multiple descriptors or multiple vectors at once, not both (or the same number of both, than it will must be a multidot)

        n = np.max((R_d_desc.shape[0], vecs.shape[0]))
        i, j = self.tril_indices

        out = np.zeros((n, self.n_atoms, self.n_atoms, 3))
        out[:, i, j, :] = R_d_desc * vecs[..., None]
        out[:, j, i, :] = -out[:, i, j, :]
        return out.sum(axis=1).reshape(n, -1)

        # if out is None or out.shape != (n, self.n_atoms*3):
        #    out = np.zeros((n, self.n_atoms*3))

        # R_d_desc_full = np.zeros((self.n_atoms, self.n_atoms, 3))
        # for a in range(n):

        #   R_d_desc_full[i, j, :] = R_d_desc * vecs[a, :, None]
        #    R_d_desc_full[j, i, :] = -R_d_desc_full[i, j, :]
        #    out[a,:] = R_d_desc_full.sum(axis=0).ravel()

        # return out

    def d_desc_from_comp(self, R_d_desc, out=None):

        if R_d_desc.ndim == 2:
            R_d_desc = R_d_desc[None, ...]

        n = R_d_desc.shape[0]
        i, j = self.tril_indices

        if out is None:
            out = np.zeros((n, self.dim, self.n_atoms, 3))
        else:
            out = out.reshape(n, self.dim, self.n_atoms, 3)

        out[:, self.dim_range, j, :] = R_d_desc
        out[:, self.dim_range, i, :] = -R_d_desc

        return out.reshape(-1, self.dim, self.dim_i)

    def d_desc_to_comp(self, R_d_desc):

        # Add singleton dimension for single inputs.
        if R_d_desc.ndim == 2:
            R_d_desc = R_d_desc[None, ...]

        n = R_d_desc.shape[0]
        n_atoms = int(R_d_desc.shape[2] / 3)

        R_d_desc = R_d_desc.reshape(n, -1, n_atoms, 3)

        ret = np.zeros((n, n_atoms, n_atoms, 3))
        ret[:, self.M, self.A, :] = R_d_desc[:, self.d_desc_mask.ravel(), self.A, :]

        # Take the upper triangle.
        i, j = self.tril_indices
        return ret[:, i, j, :]

    @staticmethod
    def perm(perm):
        n = len(perm)
        rest = np.zeros((n, n))
        rest[np.tril_indices(n, -1)] = list(range((n ** 2 - n) // 2))
        rest = rest + rest.T
        rest = rest[perm, :]
        rest = rest[:, perm]
        return rest[np.tril_indices(n, -1)].astype(int)
