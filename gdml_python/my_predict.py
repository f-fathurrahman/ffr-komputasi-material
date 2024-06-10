import sys
import logging
import os
import psutil

import numpy as np

from my_desc import Desc, _from_r


def _predict_wkr(
    train_obj, r, r_desc_d_desc, lat_and_inv, wkr_start_stop=None, chunk_size=None
):

    print("\n ***** ENTER _predict_wkr\n")

    sig, n_perms = train_obj.sig, train_obj.n_perms
    desc_func = train_obj.desc_func
    # These are from training data
    R_desc_perms = train_obj.R_desc_perms
    R_d_desc_alpha_perms = train_obj.R_d_desc_alpha_perms
    if train_obj.alphas_E_lin is not None:
        alphas_E_lin = train_obj.alphas_E_lin

    print("r is None: ", (r is None))

    print("Training: R_desc_perms.shape = ", R_desc_perms.shape)
    print("Training: R_d_desc_alpha_perms.shape = ", R_d_desc_alpha_perms.shape)

    r_desc, r_d_desc = r_desc_d_desc or desc_func.from_R(
        r, lat_and_inv, max_processes=1
    )  # no additional forking during parallelization
    print("Input: r_desc.shape = ", r_desc.shape)
    print("Input: r_d_desc.shape = ", r_d_desc.shape)

    n_train = int(R_desc_perms.shape[0] / n_perms)
    print("n_train = ", n_train)

    wkr_start, wkr_stop = (0, n_train) if wkr_start_stop is None else wkr_start_stop
    if chunk_size is None:
        chunk_size = n_train

    dim_d = desc_func.dim
    dim_i = desc_func.dim_i
    dim_c = chunk_size * n_perms

    # Pre-allocate memory.
    diff_ab_perms = np.empty((dim_c, dim_d))
    a_x2 = np.empty((dim_c,))
    mat52_base = np.empty((dim_c,))

    # avoid divisions (slower)
    sig_inv = 1.0 / sig
    mat52_base_fact = 5.0 / (3 * sig ** 3)
    diag_scale_fact = 5.0 / sig
    sqrt5 = np.sqrt(5.0)

    E_F = np.zeros((dim_d + 1,))
    F = E_F[1:]
    # First data is energy (?) the rest are forces (using linear index)

    wkr_start *= n_perms
    wkr_stop *= n_perms

    b_start = wkr_start
    print("b_start = ", b_start) 
    print("wkr_start = ", wkr_start)
    print("wkr_stop = ", wkr_stop)
    #
    # This is loop over all training data block
    #
    for b_stop in list(range(wkr_start + dim_c, wkr_stop, dim_c)) + [wkr_stop]:
        #
        print("Loop over b_stop = ", b_stop)
        #
        rj_desc_perms = R_desc_perms[b_start:b_stop, :]
        print("rj_desc_perms.shape = ", rj_desc_perms.shape)
        #
        rj_d_desc_alpha_perms = R_d_desc_alpha_perms[b_start:b_stop, :]
        print("rj_d_desc_alpha_perms.shape = ", rj_d_desc_alpha_perms.shape)
        #
        # Resize pre-allocated memory for last iteration, if chunk_size is not a divisor of the training set size.
        # Note: It's faster to process equally sized chunks.
        c_size = b_stop - b_start
        print("c_size = ", c_size)
        print("dim_c = ", dim_c)
        if c_size < dim_c:
            diff_ab_perms = diff_ab_perms[:c_size, :]
            a_x2 = a_x2[:c_size]
            mat52_base = mat52_base[:c_size]

        print("r_desc.shape = ", r_desc.shape)
        print("rj_desc_perms.shape = ", rj_desc_perms.shape)
        print("diff_ab_perms.shape = ", diff_ab_perms.shape)
        #
        np.subtract(
            np.broadcast_to(r_desc, rj_desc_perms.shape),
            rj_desc_perms,
            out=diff_ab_perms,
        )
        norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1)

        np.exp(-norm_ab_perms * sig_inv, out=mat52_base)
        mat52_base *= mat52_base_fact
        print("rj_d_desc_alpha_perms.shape = ", rj_d_desc_alpha_perms.shape)
        print("a_x2.shape = ", a_x2.shape)
        np.einsum(
            'ji,ji->j', diff_ab_perms, rj_d_desc_alpha_perms, out=a_x2
        )  # colum wise dot product

        F += (a_x2 * mat52_base).dot(diff_ab_perms) * diag_scale_fact
        mat52_base *= norm_ab_perms + sig
        F -= mat52_base.dot(rj_d_desc_alpha_perms)

        # Note: Energies are automatically predicted with a flipped sign here (because -E are trained, instead of E)
        E_F[0] += a_x2.dot(mat52_base)

        # Note: Energies are automatically predicted with a flipped sign here (because -E are trained, instead of E)
        if train_obj.alphas_E_lin is not None:

            K_fe = diff_ab_perms * mat52_base[:, None]
            F += alphas_E_lin[b_start:b_stop].dot(K_fe)

            K_ee = (
                1 + (norm_ab_perms * sig_inv) * (1 + norm_ab_perms / (3 * sig))
            ) * np.exp(-norm_ab_perms * sig_inv)
            E_F[0] += K_ee.dot(alphas_E_lin[b_start:b_stop])

        b_start = b_stop

    out = E_F[: dim_i + 1]

    # Descriptor has less entries than 3N, need to extend size of the 'E_F' array.
    if dim_d < dim_i:
        out = np.empty((dim_i + 1,))
        out[0] = E_F[0]

    out[1:] = desc_func.vec_dot_d_desc(
        r_d_desc,
        F,
    )  # 'r_d_desc.T.dot(F)' for our special representation of 'r_d_desc'


    print("\n ***** EXIT _predict_wkr\n")

    return out




class GDMLPredict(object):
    def __init__(
        self,
        model,
        batch_size=None,
        max_memory=None,
        max_processes=None,
        log_level=None,
    ):

        self.log = logging.getLogger(__name__)
        if log_level is not None:
            self.log.setLevel(log_level)

        total_memory = psutil.virtual_memory().total // 2 ** 30  # bytes to GB)
        self.max_memory = (
            min(max_memory, total_memory) if max_memory is not None else total_memory
        )

        total_cpus = 1
        self.max_processes = (
            min(max_processes, total_cpus) if max_processes is not None else total_cpus
        )

        if 'type' not in model or not (model['type'] == 'm' or model['type'] == b'm'):
            self.log.critical('The provided data structure is not a valid model.')
            sys.exit()

        self.n_atoms = model['z'].shape[0]

        self.desc_func = Desc(self.n_atoms, max_processes=max_processes)

        # Cache for iterative training mode.
        self.R_desc = None
        self.R_d_desc = None

        self.lat_and_inv = (
            (model['lattice'], np.linalg.inv(model['lattice']))
            if 'lattice' in model
            else None
        )

        self.n_train = model['R_desc'].shape[1]
        self.sig = model["sig"]

        self.std = model['std'] if 'std' in model else 1.0
        self.c = model['c']

        self.n_perms = model['perms'].shape[0]

        self.tril_perms_lin = model['tril_perms_lin']

        # Precompute permuted training descriptors and its first derivatives multiplied with the coefficients.

        # ffr: make R_desc_perms and R_d_desc_alpha_perms properties 
        self.R_desc_perms = (
            np.tile(model['R_desc'].T, self.n_perms)[:, self.tril_perms_lin]
            .reshape(self.n_train, self.n_perms, -1, order='F')
            .reshape(self.n_train * self.n_perms, -1)
        )

        self.R_d_desc_alpha_perms = (
            np.tile(model['R_d_desc_alpha'], self.n_perms)[:, self.tril_perms_lin]
            .reshape(self.n_train, self.n_perms, -1, order='F')
            .reshape(self.n_train * self.n_perms, -1)
        )

        if 'alphas_E' in model:
            self.alphas_E_lin = np.tile(model['alphas_E'][:, None], (1, self.n_perms)).ravel()
        else:
            self.alphas_E_lin = None

        self._set_chunk_size(batch_size)



    def __del__(self):
        pass


    def set_R_desc(self, R_desc):
        self.R_desc = R_desc

    def set_R_d_desc(self, R_d_desc):
        self.R_d_desc = R_d_desc


    def set_alphas(self, alphas_F, alphas_E=None):

        assert self.R_d_desc is not None

        dim_i = self.desc.dim_i
        R_d_desc_alpha = self.desc.d_desc_dot_vec(
            self.R_d_desc, alphas_F.reshape(-1, dim_i)
        )

        R_d_desc_alpha_perms_new = np.tile(R_d_desc_alpha, self.n_perms)[
            :, self.tril_perms_lin
        ].reshape(self.n_train, self.n_perms, -1, order='F')

        #R_d_desc_alpha_perms = R_d_desc_alpha_perms
        np.copyto(self.R_d_desc_alpha_perms, R_d_desc_alpha_perms_new.ravel())

        if alphas_E is not None:

            alphas_E_lin_new = np.tile(
                alphas_E[:, None], (1, self.n_perms)
            ).ravel()

            alphas_E_lin = np.frombuffer(alphas_E_lin)
            np.copyto(alphas_E_lin, alphas_E_lin_new)



    def _set_num_workers(
        self, num_workers=None, force_reset=False
    ):
        if force_reset or self.num_workers is not num_workers:
            self.num_workers = 0

        # Data ranges for processes
        if self.bulk_mp or self.num_workers < 2:
            # wkr_starts = [self.n_train]
            wkr_starts = [0]
        else:
            wkr_starts = list(
                range(
                    0,
                    self.n_train,
                    int(np.ceil(float(self.n_train) / self.num_workers)),
                )
            )
        wkr_stops = wkr_starts[1:] + [self.n_train]

        self.wkr_starts_stops = list(zip(wkr_starts, wkr_stops))

    def _set_chunk_size(self, chunk_size=None):

        if chunk_size is None:
            chunk_size = self.n_train

        self.chunk_size = chunk_size



    def _set_bulk_mp(self, bulk_mp=False):
        bulk_mp = bool(bulk_mp)
        if self.bulk_mp is not bulk_mp:
            self.bulk_mp = bulk_mp
            # Reset data ranges for processes stored in 'wkr_starts_stops'
            self._set_num_workers(self.num_workers)

    def get_GPU_batch(self):
        pass

    


    def predict(self, R=None, return_E=True):

        print("---------------------------")
        print("ENTER GDMLPredict.predict()")
        print("---------------------------")

        # Add singleton dimension if input is (,3N).
        if R is not None and R.ndim == 1:
            R = R[None, :]

        # Use precomputed descriptors in training mode.
        is_desc_in_cache = self.R_desc is not None and self.R_d_desc is not None
        print("is_desc_in_cache = ", is_desc_in_cache)

        if R is None and not is_desc_in_cache:
            self.log.critical(
                'A reference to the training geometry descriptors and Jacobians needs to be set for this function to work without arguments.'
            )
            print()
            os._exit(1)

        assert is_desc_in_cache or R is not None

        dim_i = 3 * self.n_atoms
        n_pred = self.R_desc.shape[0] if R is None else R.shape[0]

        E_F = np.empty((n_pred, dim_i + 1))

        print("n_pred = ", n_pred)
        for i in range(n_pred):

            if R is not None:
                print("R[i] = ", R[i])

            if is_desc_in_cache:
                # In case of using training data this will always be the true
                print("desc is in cache")
                r_desc, r_d_desc = self.R_desc[i], self.R_d_desc[i]
            else:
                print("desc is not is cache, calculate descriptor and its Jacobian")
                #r_desc, r_d_desc = self.desc_func.from_R(R[i], self.lat_and_inv)
                r_desc, r_d_desc = _from_r(R[i])

            E_F[i, :] = _predict_wkr(
                self, None, (r_desc, r_d_desc),
                self.lat_and_inv,
                chunk_size=self.chunk_size)

        print("self.std = ", self.std)
        print("self.c = ", self.c)

        E_F *= self.std
        F = E_F[:, 1:]
        E = E_F[:, 0] + self.c

        ret = (F,)
        if return_E:
            ret = (E,) + ret

        print("---------------------------")
        print("EXIT GDMLPredict.predict()")
        print("---------------------------")

        return ret
