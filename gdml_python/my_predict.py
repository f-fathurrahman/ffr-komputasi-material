import sys
import logging
import os
import psutil

import timeit

import numpy as np

from my_desc import Desc, _from_r


def _predict_wkr(
    train_obj, r, r_desc_d_desc, lat_and_inv, wkr_start_stop=None, chunk_size=None
):
    """
    Compute (part) of a prediction.

    Every prediction is a linear combination involving the training points used for
    this model. This function evalutates that combination for the range specified by
    `wkr_start_stop`. This workload can optionally be processed in chunks,
    which can be faster as it requires less memory to be allocated.

    Note
    ----
        It is sufficient to provide either the parameter `r` or `r_desc_d_desc`.
        The other one can be set to `None`.

    Parameters
    ----------
            r : :obj:`numpy.ndarray`
                    An array of size 3N containing the Cartesian
                    coordinates of each atom in the molecule.
            r_desc_d_desc : tuple of :obj:`numpy.ndarray`
                    A tuple made up of:
                        (1) An array of size D containing the descriptors
                        of dimension D for the molecule.
                        (2) An array of size D x 3N containing the
                        descriptor Jacobian for the molecules. It has dimension
                        D with 3N partial derivatives with respect to the 3N
                        Cartesian coordinates of each atom.
            lat_and_inv : tuple of :obj:`numpy.ndarray`
                    Tuple of 3 x 3 matrix containing lattice vectors as columns and
                    its inverse.
            glob_id : int
                    Identifier of the global namespace that this
                    function is supposed to be using (zero if only one
                    instance of this class exists at the same time).
            wkr_start_stop : tuple of int, optional
                    Range defined by the indices of first and last (exclusive)
                    sum element. The full prediction is generated if this parameter
                    is not specified.
            chunk_size : int, optional
                    Chunk size. The whole linear combination is evaluated in a large
                    vector operation instead of looping over smaller chunks if this
                    parameter is left unspecified.

    Returns
    -------
            :obj:`numpy.ndarray`
                    Partial prediction of all force components and
                    energy (appended to array as last element).
    """

    sig, n_perms = train_obj.sig, train_obj.n_perms

    desc_func = train_obj.desc_func

    R_desc_perms = train_obj.R_desc_perms
    R_d_desc_alpha_perms = train_obj.R_d_desc_alpha_perms

    if train_obj.alphas_E_lin is not None:
        alphas_E_lin = train_obj.alphas_E_lin

    r_desc, r_d_desc = r_desc_d_desc or desc_func.from_R(
        r, lat_and_inv, max_processes=1
    )  # no additional forking during parallelization

    n_train = int(R_desc_perms.shape[0] / n_perms)

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

    wkr_start *= n_perms
    wkr_stop *= n_perms

    b_start = wkr_start
    for b_stop in list(range(wkr_start + dim_c, wkr_stop, dim_c)) + [wkr_stop]:

        rj_desc_perms = R_desc_perms[b_start:b_stop, :]
        rj_d_desc_alpha_perms = R_d_desc_alpha_perms[b_start:b_stop, :]

        # Resize pre-allocated memory for last iteration, if chunk_size is not a divisor of the training set size.
        # Note: It's faster to process equally sized chunks.
        c_size = b_stop - b_start
        if c_size < dim_c:
            diff_ab_perms = diff_ab_perms[:c_size, :]
            a_x2 = a_x2[:c_size]
            mat52_base = mat52_base[:c_size]

        print(r_desc.shape)
        print("r_desc = ", r_desc)
        print(rj_desc_perms.shape)
        print(diff_ab_perms.shape)
        np.subtract(
            np.broadcast_to(r_desc, rj_desc_perms.shape),
            rj_desc_perms,
            out=diff_ab_perms,
        )
        print("Pass here 126")
        norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1)

        np.exp(-norm_ab_perms * sig_inv, out=mat52_base)
        mat52_base *= mat52_base_fact
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


    def _set_batch_size(self, batch_size=None):  # deprecated
        self._set_chunk_size(batch_size)



    def _set_bulk_mp(self, bulk_mp=False):
        bulk_mp = bool(bulk_mp)
        if self.bulk_mp is not bulk_mp:
            self.bulk_mp = bulk_mp
            # Reset data ranges for processes stored in 'wkr_starts_stops'
            self._set_num_workers(self.num_workers)


    def set_opt_num_workers_and_batch_size_fast(self, n_bulk=1, n_reps=1):  # deprecated
        self.prepare_parallel(n_bulk, n_reps)



    def prepare_parallel(
        self, n_bulk=1, n_reps=1, return_is_from_cache=False
    ):

        # Retrieve cached benchmark results, if available.
        bmark_result = self._load_cached_bmark_result(n_bulk)
        if bmark_result is not None:

            num_workers, chunk_size, bulk_mp, gps = bmark_result

            self._set_chunk_size(chunk_size)
            self._set_num_workers(num_workers)
            self._set_bulk_mp(bulk_mp)

            if return_is_from_cache:
                is_from_cache = True
                return gps, is_from_cache
            else:
                return gps

        warm_up_done = False

        best_results = []
        last_i = None

        best_gps = 0
        gps_min = 0.0

        best_params = None

        r_dummy = np.random.rand(n_bulk, self.n_atoms * 3)

        def _dummy_predict():
            self.predict(r_dummy)

        bulk_mp_rng = [True, False] if n_bulk > 1 else [False]
        for bulk_mp in bulk_mp_rng:
            self._set_bulk_mp(bulk_mp)

            if bulk_mp is False:
                last_i = 0

            num_workers_rng = list(range(0, self.max_processes))
            if bulk_mp:
                num_workers_rng.reverse()  # benchmark converges faster this way

            # num_workers_rng_sizes = [batch_size for batch_size in batch_size_rng if min_batch_size % batch_size == 0]

            # for num_workers in range(min_num_workers,self.max_processes+1):
            for num_workers in num_workers_rng:
                if not bulk_mp and num_workers != 0 and self.n_train % num_workers != 0:
                    continue

                self._set_num_workers(num_workers)

                best_gps = 0
                gps_rng = (np.inf, 0.0)  # min and max per num_workers

                min_chunk_size = (
                    min(self.n_train, n_bulk)
                    if bulk_mp or num_workers < 2
                    else int(np.ceil(self.n_train / num_workers))
                )
                chunk_size_rng = list(range(min_chunk_size, 0, -1))

                chunk_size_rng_sizes = [
                    chunk_size
                    for chunk_size in chunk_size_rng
                    if min_chunk_size % chunk_size == 0
                ]

                # print('batch_size_rng_sizes ' + str(bulk_mp))
                # print(batch_size_rng_sizes)

                i_done = 0
                i_dir = 1
                i = 0 if last_i is None else last_i
                # i = 0

                # print(batch_size_rng_sizes)
                while i >= 0 and i < len(chunk_size_rng_sizes):

                    chunk_size = chunk_size_rng_sizes[i]
                    self._set_chunk_size(chunk_size)

                    i_done += 1

                    if warm_up_done == False:
                        timeit.timeit(_dummy_predict, number=10)
                        warm_up_done = True

                    gps = n_bulk * n_reps / timeit.timeit(_dummy_predict, number=n_reps)

                    # print(
                    #  '{:2d}@{:d} {:d} | {:7.2f} gps'.format(
                    #      num_workers, chunk_size, bulk_mp, gps
                    #  )
                    # )

                    gps_rng = (
                        min(gps_rng[0], gps),
                        max(gps_rng[1], gps),
                    )  # min and max per num_workers

                    # gps still going up?
                    # AND: gps not lower than the lowest overall?
                    # if gps < best_gps and gps >= gps_min:
                    if gps < best_gps:
                        if (
                            i_dir > 0
                            and i_done == 2
                            and chunk_size
                            != chunk_size_rng_sizes[
                                1
                            ]  # there is no point in turning if this is the second batch size in the range
                        ):  # do we turn?
                            i -= 2 * i_dir
                            i_dir = -1
                            # print('><')
                            continue
                        else:
                            if chunk_size == chunk_size_rng_sizes[1]:
                                i -= 1 * i_dir
                            # print('>>break ' + str(i_done))
                            break
                    else:
                        best_gps = gps
                        best_params = num_workers, chunk_size, bulk_mp

                    if (
                        not bulk_mp and n_bulk > 1
                    ):  # stop search early when multiple cpus are available and the 1 cpu case is tested
                        if (
                            gps < gps_min
                        ):  # if the batch size run is lower than the lowest overall, stop right here
                            # print('breaking here')
                            break

                    i += 1 * i_dir

                last_i = i - 1 * i_dir
                i_dir = 1

                if len(best_results) > 0:
                    overall_best_gps = max(best_results, key=lambda x: x[1])[1]
                    if best_gps < overall_best_gps:
                        # print('breaking, because best of last test was worse than overall best so far')
                        break

                    # if best_gps < gps_min:
                    #    print('breaking here3')
                    #    break

                gps_min = gps_rng[0]  # FIX me: is this the overall min?
                # print ('gps_min ' + str(gps_min))

                # print ('best_gps')
                # print (best_gps)

                best_results.append(
                    (best_params, best_gps)
                )  # best results per num_workers

        (num_workers, chunk_size, bulk_mp), gps = max(best_results, key=lambda x: x[1])

        # Cache benchmark results.
        self._save_cached_bmark_result(n_bulk, num_workers, chunk_size, bulk_mp, gps)

        self._set_chunk_size(chunk_size)
        self._set_num_workers(num_workers)
        self._set_bulk_mp(bulk_mp)

        if return_is_from_cache:
            is_from_cache = False
            return gps, is_from_cache
        else:
            return gps

    def _save_cached_bmark_result(self, n_bulk, num_workers, chunk_size, bulk_mp, gps):

        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        bmark_file = '_bmark_cache.npz'
        bmark_path = os.path.join(pkg_dir, bmark_file)

        bkey = '{}-{}-{}-{}'.format(
            self.n_atoms, self.n_train, n_bulk, self.max_processes
        )

        if os.path.exists(bmark_path):

            with np.load(bmark_path, allow_pickle=True) as bmark:
                bmark = dict(bmark)

                bmark['runs'] = np.append(bmark['runs'], bkey)
                bmark['num_workers'] = np.append(bmark['num_workers'], num_workers)
                bmark['batch_size'] = np.append(bmark['batch_size'], chunk_size)
                bmark['bulk_mp'] = np.append(bmark['bulk_mp'], bulk_mp)
                bmark['gps'] = np.append(bmark['gps'], gps)
        else:
            bmark = {
                'code_version': '0.0.1.ffr',
                'runs': [bkey],
                'gps': [gps],
                'num_workers': [num_workers],
                'batch_size': [chunk_size],
                'bulk_mp': [bulk_mp],
            }

        np.savez_compressed(bmark_path, **bmark)

    def _load_cached_bmark_result(self, n_bulk):

        pkg_dir = os.path.dirname(os.path.abspath(__file__))
        bmark_file = '_bmark_cache.npz'
        bmark_path = os.path.join(pkg_dir, bmark_file)

        bkey = '{}-{}-{}-{}'.format(
            self.n_atoms, self.n_train, n_bulk, self.max_processes
        )

        if not os.path.exists(bmark_path):
            return None

        with np.load(bmark_path, allow_pickle=True) as bmark:

            # Keep collecting benchmark runs, until we have at least three.
            run_idxs = np.where(bmark['runs'] == bkey)[0]
            if len(run_idxs) >= 3:

                config_keys = []
                for run_idx in run_idxs:
                    config_keys.append(
                        '{}-{}-{}'.format(
                            bmark['num_workers'][run_idx],
                            bmark['batch_size'][run_idx],
                            bmark['bulk_mp'][run_idx],
                        )
                    )

                values, uinverse = np.unique(config_keys, return_index=True)

                best_mean = -1
                best_gps = 0
                for i, config_key in enumerate(zip(values, uinverse)):
                    mean_gps = np.mean(
                        bmark['gps'][
                            np.where(np.array(config_keys) == config_key[0])[0]
                        ]
                    )

                    if best_gps == 0 or best_gps < mean_gps:
                        best_mean = i
                        best_gps = mean_gps

                best_idx = run_idxs[uinverse[best_mean]]
                num_workers = bmark['num_workers'][best_idx]
                chunk_size = bmark['batch_size'][best_idx]
                bulk_mp = bmark['bulk_mp'][best_idx]

                return num_workers, chunk_size, bulk_mp, best_gps

        return None

    def get_GPU_batch(self):
        pass

    def predict(self, R=None, return_E=True):
        # Add singleton dimension if input is (,3N).
        if R is not None and R.ndim == 1:
            R = R[None, :]

        # Use precomputed descriptors in training mode.
        is_desc_in_cache = self.R_desc is not None and self.R_d_desc is not None

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

            print("R[i] = ", R[i])

            if is_desc_in_cache:
                r_desc, r_d_desc = self.R_desc[i], self.R_d_desc[i]
            else:
                #r_desc, r_d_desc = self.desc_func.from_R(R[i], self.lat_and_inv)
                r_desc, r_d_desc = _from_r(R[i])

            E_F[i, :] = _predict_wkr(
                self, None, (r_desc, r_d_desc),
                self.lat_and_inv,
                chunk_size=self.chunk_size)

        E_F *= self.std
        F = E_F[:, 1:]
        E = E_F[:, 0] + self.c

        ret = (F,)
        if return_E:
            ret = (E,) + ret

        return ret
