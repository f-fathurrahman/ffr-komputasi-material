import os
import logging
import psutil

import timeit
from functools import partial

import numpy as np

from . import __version__, DONE, NOT_DONE

from my_predict import GDMLPredict
from my_desc import Desc


import hashlib

def io_dataset_md5(dataset):

    md5_hash = hashlib.md5()

    keys = ['z', 'R']
    if 'E' in dataset:
        keys.append('E')
    keys.append('F')

    for k in keys:
        d = dataset[k]
        if type(d) is np.ndarray:
            d = d.ravel()
        md5_hash.update(hashlib.md5(d).digest())

    return md5_hash.hexdigest().encode('utf-8')





def _assemble_kernel_mat_wkr(
    K, desc_func, R_desc, R_d_desc,
    j, tril_perms_lin, sig, use_E_cstr=False, exploit_sym=False, cols_m_limit=None
):

    n_train, dim_d = R_d_desc.shape[:2]
    n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)
    dim_i = 3 * n_atoms
    n_perms = int(len(tril_perms_lin) / dim_d)

    if type(j) is tuple:  # Selective/"fancy" indexing
        (
            K_j,
            j,
            keep_idxs_3n,
        ) = j  # (block index in final K, block index global, indices of partials within block)
        blk_j = slice(K_j, K_j + len(keep_idxs_3n))

    else:  # Sequential indexing
        blk_j = slice(j * dim_i, (j + 1) * dim_i)
        keep_idxs_3n = slice(None)  # same as [:]

    # TODO: document this exception
    if use_E_cstr and not (cols_m_limit is None or cols_m_limit == n_train):
        raise ValueError(
            '\'use_E_cstr\'- and \'cols_m_limit\'-parameters are mutually exclusive!'
        )

    # Create permutated variants of 'rj_desc' and 'rj_d_desc'.
    rj_desc_perms = np.reshape(
        np.tile(R_desc[j, :], n_perms)[tril_perms_lin], (n_perms, -1), order='F'
    )

    rj_d_desc = desc_func.d_desc_from_comp(R_d_desc[j, :, :])[0][
        :, keep_idxs_3n
    ]  # convert descriptor back to full representation

    rj_d_desc_perms = np.reshape(
        np.tile(rj_d_desc.T, n_perms)[:, tril_perms_lin], (-1, dim_d, n_perms)
    )

    mat52_base_div = 3 * sig ** 4
    sqrt5 = np.sqrt(5.0)
    sig_pow2 = sig ** 2

    dim_i_keep = rj_d_desc.shape[1]
    diff_ab_outer_perms = np.empty((dim_d, dim_i_keep))
    diff_ab_perms = np.empty((n_perms, dim_d))
    ri_d_desc = np.zeros((1, dim_d, dim_i))  # must be zeros!
    k = np.empty((dim_i, dim_i_keep))

    for i in range(j if exploit_sym else 0, n_train):

        blk_i = slice(i * dim_i, (i + 1) * dim_i)

        # diff_ab_perms = R_desc[i, :] - rj_desc_perms
        np.subtract(R_desc[i, :], rj_desc_perms, out=diff_ab_perms)

        norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1)
        mat52_base_perms = np.exp(-norm_ab_perms / sig) / mat52_base_div * 5

        # diff_ab_outer_perms = 5 * np.einsum(
        #    'ki,kj->ij',
        #    diff_ab_perms * mat52_base_perms[:, None],
        #    np.einsum('ik,jki -> ij', diff_ab_perms, rj_d_desc_perms)
        # )
        np.einsum(
            'ki,kj->ij',
            diff_ab_perms * mat52_base_perms[:, None] * 5,
            np.einsum('ki,jik -> kj', diff_ab_perms, rj_d_desc_perms),
            out=diff_ab_outer_perms,
        )

        diff_ab_outer_perms -= np.einsum(
            'ikj,j->ki',
            rj_d_desc_perms,
            (sig_pow2 + sig * norm_ab_perms) * mat52_base_perms,
        )

        # ri_d_desc = desc_func.d_desc_from_comp(R_d_desc[i, :, :])[0]
        desc_func.d_desc_from_comp(R_d_desc[i, :, :], out=ri_d_desc)

        # K[blk_i, blk_j] = ri_d_desc[0].T.dot(diff_ab_outer_perms)
        np.dot(ri_d_desc[0].T, diff_ab_outer_perms, out=k)
        K[blk_i, blk_j] = k

        # print(k - k2)

        if exploit_sym and (
            cols_m_limit is None or i < cols_m_limit
        ):  # this will never be called with 'keep_idxs_3n' set to anything else than [:]
            K[blk_j, blk_i] = K[blk_i, blk_j].T

    if use_E_cstr:

        E_off = K.shape[0] - n_train, K.shape[1] - n_train
        blk_j_full = slice(j * dim_i, (j + 1) * dim_i)
        for i in range(n_train):

            diff_ab_perms = R_desc[i, :] - rj_desc_perms
            norm_ab_perms = sqrt5 * np.linalg.norm(diff_ab_perms, axis=1)

            K_fe = (
                5
                * diff_ab_perms
                / (3 * sig ** 3)
                * (norm_ab_perms[:, None] + sig)
                * np.exp(-norm_ab_perms / sig)[:, None]
            )
            K_fe = -np.einsum('ik,jki -> j', K_fe, rj_d_desc_perms)
            K[blk_j_full, E_off[1] + i] = K_fe  # vertical
            K[E_off[0] + i, blk_j] = K_fe[keep_idxs_3n]  # lower horizontal

            K[E_off[0] + i, E_off[1] + j] = K[E_off[0] + j, E_off[1] + i] = -(
                1 + (norm_ab_perms / sig) * (1 + norm_ab_perms / (3 * sig))
            ).dot(np.exp(-norm_ab_perms / sig))

    return blk_j.stop - blk_j.start




class GDMLTrain(object):
    def __init__(self, max_memory=None, max_processes=None, use_torch=False):
        self.log = logging.getLogger(__name__)
        total_memory = psutil.virtual_memory().total // 2 ** 30  # bytes to GB)
        self._max_memory = (
            min(max_memory, total_memory) if max_memory is not None else total_memory
        )
        self._max_processes = 1


    def __del__(self):
        pass



    def create_task(
        self,
        train_dataset,
        n_train,
        valid_dataset,
        n_valid,
        sig,
        lam=1e-10,
        perms=None,
        use_sym=True,
        use_E=True,
        use_E_cstr=False,
        callback=None,  # TODO: document me
    ):

        if use_E and 'E' not in train_dataset:
            raise ValueError(
                'No energy labels found in dataset!\n'
                + 'By default, force fields are always reconstructed including the\n'
                + 'corresponding potential energy surface (this can be turned off).\n'
                + 'However, the energy labels are missing in the provided dataset.\n'
            )

        use_E_cstr = use_E and use_E_cstr
        print("use_E_cstr = ", use_E_cstr)

        n_atoms = train_dataset['R'].shape[1]
        print("n_atoms = ", n_atoms)

        if callback is not None:
            callback = partial(callback, disp_str='Hashing dataset(s)')
            callback(NOT_DONE)

        md5_train = io_dataset_md5(train_dataset)
        md5_valid = io_dataset_md5(valid_dataset)

        if callback is not None:
            callback(DONE)

        if callback is not None:
            callback = partial(
                callback, disp_str='Sampling training and validation subsets'
            )
            callback(NOT_DONE)

        if 'E' in train_dataset:
            idxs_train = self.draw_strat_sample(train_dataset['E'], n_train)
        else:
            idxs_train = np.random.choice(
                np.arange(train_dataset['F'].shape[0]),
                n_train,
                replace=False,
            )

        excl_idxs = (
            idxs_train if md5_train == md5_valid else np.array([], dtype=np.uint)
        )

        if 'E' in valid_dataset:
            idxs_valid = self.draw_strat_sample(
                valid_dataset['E'],
                n_valid,
                excl_idxs=excl_idxs,
            )
        else:
            idxs_valid_cands = np.setdiff1d(
                np.arange(valid_dataset['F'].shape[0]), excl_idxs, assume_unique=True
            )
            idxs_valid = np.random.choice(idxs_valid_cands, n_valid, replace=False)

        if callback is not None:
            callback(DONE)

        R_train = train_dataset['R'][idxs_train, :, :]
        task = {
            'type': 't',
            'code_version': __version__,
            'dataset_name': train_dataset['name'].astype(str),
            'dataset_theory': train_dataset['theory'].astype(str),
            'z': train_dataset['z'],
            'R_train': R_train,
            'F_train': train_dataset['F'][idxs_train, :, :],
            'idxs_train': idxs_train,
            'md5_train': md5_train,
            'idxs_valid': idxs_valid,
            'md5_valid': md5_valid,
            'sig': sig,
            'lam': lam,
            'use_E': use_E,
            'use_E_cstr': use_E_cstr,
            'use_sym': use_sym,
        }

        if use_E:
            task['E_train'] = train_dataset['E'][idxs_train]

        lat_and_inv = None
        if 'lattice' in train_dataset:
            task['lattice'] = train_dataset['lattice']

            try:
                lat_and_inv = (task['lattice'], np.linalg.inv(task['lattice']))
            except np.linalg.LinAlgError:
                raise ValueError(  # TODO: Document me
                    'Provided dataset contains invalid lattice vectors (not invertible). Note: Only rank 3 lattice vector matrices are supported.'
                )

        if 'r_unit' in train_dataset and 'e_unit' in train_dataset:
            task['r_unit'] = train_dataset['r_unit']
            task['e_unit'] = train_dataset['e_unit']

        if use_sym:
            print("use_symm is disabled")

        else:
            task['perms'] = np.arange(train_dataset['R'].shape[1])[
                None, :
            ]  # no symmetries

        return task




    def create_task_from_model(self, model, dataset):

        idxs_train = model['idxs_train']
        R_train = dataset['R'][idxs_train, :, :]
        F_train = dataset['F'][idxs_train, :, :]

        use_E = 'e_err' in model
        use_E_cstr = 'alphas_E' in model
        use_sym = model['perms'].shape[0] > 1

        task = {
            'type': 't',
            'code_version': __version__,
            'dataset_name': model['dataset_name'],
            'dataset_theory': model['dataset_theory'],
            'z': model['z'],
            'R_train': R_train,
            'F_train': F_train,
            'idxs_train': idxs_train,
            'md5_train': model['md5_train'],
            'idxs_valid': model['idxs_valid'],
            'md5_valid': model['md5_valid'],
            'sig': model['sig'],
            'lam': model['lam'],
            'use_E': model['use_E'],
            'use_E_cstr': use_E_cstr,
            'use_sym': use_sym,
            'perms': model['perms'],
        }

        if use_E:
            task['E_train'] = dataset['E'][idxs_train]

        if 'lattice' in model:
            task['lattice'] = model['lattice']

        if 'r_unit' in model and 'e_unit' in model:
            task['r_unit'] = model['r_unit']
            task['e_unit'] = model['e_unit']

        if 'alphas_F' in model:
            task['alphas0_F'] = model['alphas_F']

        if 'alphas_E' in model:
            task['alphas0_E'] = model['alphas_E']

        if 'solver_iters' in model:
            task['solver_iters'] = model['solver_iters']

        if 'inducing_pts_idxs' in model:
            task['inducing_pts_idxs'] = model['inducing_pts_idxs']

        return task

    def create_model(
        self,
        task,
        solver,
        R_desc,
        R_d_desc,
        tril_perms_lin,
        std,
        alphas_F,
        alphas_E=None,
    ):
        
        n_train, dim_d = R_d_desc.shape[:2]
        n_atoms = int((1 + np.sqrt(8 * dim_d + 1)) / 2)

        desc = Desc(
            n_atoms,
            max_processes=self._max_processes,
        )

        dim_i = desc.dim_i
        R_d_desc_alpha = desc.d_desc_dot_vec(R_d_desc, alphas_F.reshape(-1, dim_i))

        model = {
            'type': 'm',
            'code_version': __version__,
            'dataset_name': task['dataset_name'],
            'dataset_theory': task['dataset_theory'],
            'solver_name': solver,
            'z': task['z'],
            'idxs_train': task['idxs_train'],
            'md5_train': task['md5_train'],
            'idxs_valid': task['idxs_valid'],
            'md5_valid': task['md5_valid'],
            'n_test': 0,
            'md5_test': None,
            'f_err': {'mae': np.nan, 'rmse': np.nan},
            'R_desc': R_desc.T,
            'R_d_desc_alpha': R_d_desc_alpha,
            'c': 0.0,
            'std': std,
            'sig': task['sig'],
            'lam': task['lam'],
            'alphas_F': alphas_F,
            'perms': task['perms'],
            'tril_perms_lin': tril_perms_lin,
            'use_E': task['use_E'],
        }

        if task['use_E']:
            model['e_err'] = {'mae': np.nan, 'rmse': np.nan}

            if task['use_E_cstr']:
                model['alphas_E'] = alphas_E

        if 'lattice' in task:
            model['lattice'] = task['lattice']

        if 'r_unit' in task and 'e_unit' in task:
            model['r_unit'] = task['r_unit']
            model['e_unit'] = task['e_unit']

        return model


    def train(
        self,
        task,
        save_progr_callback=None,  # TODO: document me
        callback=None,
    ):

        task = dict(task)  # make mutable

        n_train, n_atoms = task['R_train'].shape[:2]

        desc = Desc(
            n_atoms,
            max_processes=self._max_processes,
        )

        n_perms = task['perms'].shape[0]
        tril_perms = np.array([Desc.perm(p) for p in task['perms']])

        dim_i = 3 * n_atoms
        dim_d = desc.dim

        perm_offsets = np.arange(n_perms)[:, None] * dim_d
        tril_perms_lin = (tril_perms + perm_offsets).flatten('F')

        # TODO: check if all atoms are in span of lattice vectors, otherwise suggest that
        # rows and columns might have been switched.
        lat_and_inv = None
        if 'lattice' in task:
            try:
                lat_and_inv = (task['lattice'], np.linalg.inv(task['lattice']))
            except np.linalg.LinAlgError:
                raise ValueError(  # TODO: Document me
                    'Provided dataset contains invalid lattice vectors (not invertible). Note: Only rank 3 lattice vector matrices are supported.'
                )

            # # TODO: check if all atoms are within unit cell
            # for r in task['R_train']:
            #    r_lat = lat_and_inv[1].dot(r.T)
            #    if not (r_lat >= 0).all():
            #         raise ValueError( # TODO: Document me
            #            'Some atoms appear outside of the unit cell! Please check lattice vectors in dataset file.'
            #         )
            #        #pass

        R = task['R_train'].reshape(n_train, -1)
        R_desc, R_d_desc = desc.from_R(
            R,
            lat_and_inv=lat_and_inv,
            callback=partial(
                callback, disp_str='Generating descriptors and their Jacobians'
            )
            if callback is not None
            else None,
        )

        # Generate label vector.
        E_train_mean = None
        y = task['F_train'].ravel().copy()
        if task['use_E'] and task['use_E_cstr']:
            E_train = task['E_train'].ravel().copy()
            E_train_mean = np.mean(E_train)

            y = np.hstack((y, -E_train + E_train_mean))
            # y = np.hstack((n*Ft, (1-n)*Et))
        y_std = np.std(y)
        print("y_std = ", y_std)
        y /= y_std

        max_memory_bytes = self._max_memory * 1024 ** 3

        # Memory cost of analytic solver
        est_bytes_analytic = Analytic.est_memory_requirement(n_train, n_atoms)

        # Memory overhead (solver independent)
        est_bytes_overhead = y.nbytes
        est_bytes_overhead += R.nbytes
        est_bytes_overhead += R_desc.nbytes
        est_bytes_overhead += R_d_desc.nbytes

        solver_keys = {}

        use_analytic_solver = (
            est_bytes_analytic + est_bytes_overhead
        ) < max_memory_bytes

        print("use_analytic_solver = ", use_analytic_solver)

        # Fall back to analytic solver, if iterative solver file is missing.
        base_path = os.path.dirname(os.path.abspath(__file__))
        iter_solver_path = os.path.join(base_path, 'solvers/iterative.py')
        if not os.path.exists(iter_solver_path):
            self.log.debug('Iterative solver not installed.')
            use_analytic_solver = True

        # use_analytic_solver = False # remove me!

        if use_analytic_solver:

            self.log.info(
                'Using analytic solver (expected memory requirement: ~{})'.format(
                    ui.gen_memory_str(est_bytes_analytic + est_bytes_overhead)
                )
            )

            print("Using analytic solver")
            analytic = Analytic(self, desc, callback=callback)
            alphas = analytic.solve(task, R_desc, R_d_desc, tril_perms_lin, y)

        else:

            max_n_inducing_pts = Iterative.max_n_inducing_pts(
                n_train, n_atoms, max_memory_bytes
            )
            est_bytes_iterative = Iterative.est_memory_requirement(
                n_train, max_n_inducing_pts, n_atoms
            )

            self.log.info(
                'Using iterative solver (expected memory requirement: ~{})'.format(
                    ui.gen_memory_str(est_bytes_iterative + est_bytes_overhead)
                )
            )

            alphas_F = task['alphas0_F'] if 'alphas0_F' in task else None
            alphas_E = task['alphas0_E'] if 'alphas0_E' in task else None

            iterative = Iterative(
                self,
                desc,
                self._max_memory,
                self._max_processes,
                self._use_torch,
                callback=callback,
            )
            (
                alphas,
                solver_keys['solver_tol'],
                solver_keys[
                    'solver_iters'
                ],  # number of iterations performed (cg solver)
                solver_keys['solver_resid'],  # residual of solution
                train_rmse,
                solver_keys['inducing_pts_idxs'],
                is_conv,
            ) = iterative.solve(
                task,
                R_desc,
                R_d_desc,
                tril_perms_lin,
                y,
                y_std,
                save_progr_callback=save_progr_callback,
            )

            solver_keys['norm_y_train'] = np.linalg.norm(y)

            if not is_conv:
                print('Iterative solver did not converge!')

        alphas_E = None
        alphas_F = alphas
        if task['use_E_cstr']:
            alphas_E = alphas[-n_train:]
            alphas_F = alphas[:-n_train]


        print("After finding the parameters, create model which will be returned")
        print("solver_keys = ", solver_keys)

        model = self.create_model(
            task,
            'analytic' if use_analytic_solver else 'cg',
            R_desc,
            R_d_desc,
            tril_perms_lin,
            y_std,
            alphas_F,
            alphas_E=alphas_E,
        )
        model.update(solver_keys)

        # Recover integration constant.
        # Note: if energy constraints are included in the kernel (via 'use_E_cstr'), do not
        # compute the integration constant, but simply set it to the mean of the training energies
        # (which was subtracted from the labels before training).
        if model['use_E']:
            c = (
                self._recov_int_const(model, task, R_desc=R_desc, R_d_desc=R_d_desc)
                if E_train_mean is None
                else E_train_mean
            )
            print("Recover integration constant: c = ", c)
            if c is None:
                # Something does not seem right. Turn off energy predictions for this model, only output force predictions.
                model['use_E'] = False
            else:
                model['c'] = c

        print(">>>> EXIT GDML.train")

        return model




    def _recov_int_const(
        self, model, task, R_desc=None, R_d_desc=None
    ):  # TODO: document e_err_inconsist return

        gdml_predict = GDMLPredict(
            model,
            max_memory=self._max_memory,
            max_processes=self._max_processes,
            use_torch=self._use_torch,
            log_level=logging.CRITICAL,
        )

        gdml_predict.set_R_desc(R_desc)
        gdml_predict.set_R_d_desc(R_d_desc)

        E_pred, _ = gdml_predict.predict()
        E_ref = np.squeeze(task['E_train'])

        e_fact = np.linalg.lstsq(
            np.column_stack((E_pred, np.ones(E_ref.shape))), E_ref, rcond=-1
        )[0][0]
        corrcoef = np.corrcoef(E_ref, E_pred)[0, 1]

        # import matplotlib.pyplot as plt
        # sidx = np.argsort(E_ref)
        # plt.plot(E_ref[sidx])
        # c = np.sum(E_ref - E_pred) / E_ref.shape[0]
        # plt.plot(E_pred[sidx]+c)
        # plt.show()
        # sys.exit()

        # import matplotlib.pyplot as plt
        # sidx = np.argsort(F_ref)
        # plt.plot(F_ref[sidx])
        # c = np.sum(F_ref - F_pred) / F_ref.shape[0]
        # plt.plot(F_pred[sidx],'--')
        # plt.show()
        # sys.exit()

        if np.sign(e_fact) == -1:
            self.log.warning(
                'The provided dataset contains gradients instead of force labels (flipped sign). Please correct!\n'
                + ui.color_str('Note:', bold=True)
                + 'Note: The energy prediction accuracy of the model will thus neither be validated nor tested in the following steps!'
            )
            return None

        if corrcoef < 0.95:
            self.log.warning(
                'Inconsistent energy labels detected!\n'
                + 'The predicted energies for the training data are only weakly correlated with the reference labels (correlation coefficient {:.2f}) which indicates that the issue is most likely NOT just a unit conversion error.\n\n'.format(
                    corrcoef
                )
                + ui.color_str('Troubleshooting tips:\n', bold=True)
                + ui.wrap_indent_str(
                    '(1) ',
                    'Verify the correct correspondence between geometries and labels in the provided dataset.',
                )
                + '\n'
                + ui.wrap_indent_str(
                    '(2) ', 'Verify the consistency between energy and force labels.'
                )
                + '\n'
                + ui.wrap_indent_str('    - ', 'Correspondence correct?')
                + '\n'
                + ui.wrap_indent_str('    - ', 'Same level of theory?')
                + '\n'
                + ui.wrap_indent_str('    - ', 'Accuracy of forces?')
                + '\n'
                + ui.wrap_indent_str(
                    '(3) ',
                    'Is the training data spread too broadly (i.e. weakly sampled transitions between example clusters)?',
                )
                + '\n'
                + ui.wrap_indent_str(
                    '(4) ', 'Are there duplicate geometries in the training data?'
                )
                + '\n'
                + ui.wrap_indent_str(
                    '(5) ', 'Are there any corrupted data points (e.g. parsing errors)?'
                )
                + '\n\n'
                + ui.color_str('Note:', bold=True)
                + ' The energy prediction accuracy of the model will thus neither be validated nor tested in the following steps!'
            )
            return None

        if np.abs(e_fact - 1) > 1e-1:
            self.log.warning(
                'Different scales in energy vs. force labels detected!\n'
                + 'The integrated forces differ from the energy labels by factor ~{:.2f}, meaning that the trained model will likely fail to predict energies accurately.\n\n'.format(
                    e_fact
                )
                + ui.color_str('Troubleshooting tips:\n', bold=True)
                + ui.wrap_indent_str(
                    '(1) ', 'Verify consistency of units in energy and force labels.'
                )
                + '\n'
                + ui.wrap_indent_str(
                    '(2) ',
                    'Is the training data spread too broadly (i.e. weakly sampled transitions between example clusters)?',
                )
                + '\n\n'
                + ui.color_str('Note:', bold=True)
                + ' The energy prediction accuracy of the model will thus neither be validated nor tested in the following steps!'
            )
            return None

        # Least squares estimate for integration constant.
        return np.sum(E_ref - E_pred) / E_ref.shape[0]

    def _assemble_kernel_mat(
        self,
        R_desc,
        R_d_desc,
        tril_perms_lin,
        sig,
        desc,  # TODO: document me
        use_E_cstr=False,
        col_idxs=np.s_[:],  # TODO: document me
        alloc_extra_rows=0,  # TODO: document me
        callback=None,
    ):
        r"""
        Compute force field kernel matrix.

        The Hessian of the Matern kernel is used with n = 2 (twice
        differentiable). Each row and column consists of matrix-valued blocks,
        which encode the interaction of one training point with all others. The
        result is stored in shared memory (a global variable).

        Parameters
        ----------
            R_desc : :obj:`numpy.ndarray`
                Array containing the descriptor for each training point.
            R_d_desc : :obj:`numpy.ndarray`
                Array containing the gradient of the descriptor for
                each training point.
            tril_perms_lin : :obj:`numpy.ndarray`
                1D array containing all recovered permutations
                expanded as one large permutation to be applied to a
                tiled copy of the object to be permuted.
            sig : int
                Hyper-parameter :math:`\sigma`(kernel length scale).
            use_E_cstr : bool, optional
                True: include energy constraints in the kernel,
                False: default (s)GDML kernel.
            callback : callable, optional
                Kernel assembly progress function that takes three
                arguments:
                    current : int
                        Current progress (number of completed entries).
                    total : int
                        Task size (total number of entries to create).
                    done_str : :obj:`str`, optional
                        Once complete, this string contains the
                        time it took to assemble the kernel (seconds).
            cols_m_limit : int, optional
                Only generate the columns up to index 'cols_m_limit'. This creates
                a M*3N x cols_m_limit*3N kernel matrix, instead of M*3N x M*3N.
            cols_3n_keep_idxs : :obj:`numpy.ndarray`, optional
                Only generate columns with the given indices in the 3N x 3N
                kernel function. The resulting kernel matrix will have dimension
                M*3N x M*len(cols_3n_keep_idxs).

        Returns
        -------
            :obj:`numpy.ndarray`
                Force field kernel matrix.
        """

        global glob

        # Note: This function does not support unsorted (ascending) index arrays.
        # if not isinstance(col_idxs, slice):
        #    assert np.array_equal(col_idxs, np.sort(col_idxs))

        n_train, dim_d = R_d_desc.shape[:2]
        dim_i = 3 * int((1 + np.sqrt(8 * dim_d + 1)) / 2)

        # Determine size of kernel matrix.
        K_n_rows = n_train * dim_i
        if isinstance(col_idxs, slice):  # indexed by slice
            K_n_cols = len(range(*col_idxs.indices(K_n_rows)))
        else:  # indexed by list

            # TODO: throw exeption with description
            assert len(col_idxs) == len(set(col_idxs))  # assume no dublicate indices

            # TODO: throw exeption with description
            # Note: This function does not support unsorted (ascending) index arrays.
            assert np.array_equal(col_idxs, np.sort(col_idxs))

            K_n_cols = len(col_idxs)

        # Account for additional rows and columns due to energy constraints in the kernel matrix.
        if use_E_cstr:
            K_n_rows += n_train
            K_n_cols += n_train

        # Make sure no indices are outside of the valid range.
        if K_n_cols > K_n_rows:
            raise ValueError('Columns indexed beyond range.')

        exploit_sym = False
        cols_m_limit = None

        # Check if range is a subset of training points (as opposed to a subset of partials of multiple points).
        is_M_subset = (
            isinstance(col_idxs, slice)
            and (col_idxs.start is None or col_idxs.start % dim_i == 0)
            and (col_idxs.stop is None or col_idxs.stop % dim_i == 0)
            and col_idxs.step is None
        )
        if is_M_subset:
            M_slice_start = (
                None if col_idxs.start is None else int(col_idxs.start / dim_i)
            )
            M_slice_stop = None if col_idxs.stop is None else int(col_idxs.stop / dim_i)
            M_slice = slice(M_slice_start, M_slice_stop)

            J = range(*M_slice.indices(n_train))

            if M_slice_start is None:
                exploit_sym = True
                cols_m_limit = M_slice_stop

        else:

            if isinstance(col_idxs, slice):
                random = list(range(*col_idxs.indices(n_train * dim_i)))
            else:
                random = col_idxs

            # M - number training
            # N - number atoms

            n_idxs = np.mod(random, dim_i)
            m_idxs = (np.array(random) / dim_i).astype(int)
            m_idxs_uniq = np.unique(m_idxs)  # which points to include?

            m_n_idxs = [
                list(n_idxs[np.where(m_idxs == m_idx)]) for m_idx in m_idxs_uniq
            ]

            m_n_idxs_lens = [len(m_n_idx) for m_n_idx in m_n_idxs]

            m_n_idxs_lens.insert(0, 0)
            blk_start_idxs = list(
                np.cumsum(m_n_idxs_lens[:-1])
            )  # index within K at which each block starts

            # tupels: (block index in final K, block index global, indices of partials within block)
            J = list(zip(blk_start_idxs, m_idxs_uniq, m_n_idxs))

        if callback is not None:
            callback(0, 100)  # 0%

        if self._use_torch:
            if not _has_torch:
                raise ImportError(
                    'Optional PyTorch dependency not found! Please run \'pip install sgdml[torch]\' to install it or disable the PyTorch option.'
                )

            K = np.empty((K_n_rows + alloc_extra_rows, K_n_cols))

            if J is not list:
                J = list(J)

            global torch_assemble_done
            torch_assemble_todo, torch_assemble_done = K_n_cols, 0

            def progress_callback(done):

                global torch_assemble_done
                torch_assemble_done += done

                if callback is not None:
                    callback(
                        torch_assemble_done,
                        torch_assemble_todo,
                        newline_when_done=False,
                    )

            start = timeit.default_timer()

            torch_device = 'cuda' if torch.cuda.is_available() else 'cpu'

            R_desc_torch = torch.from_numpy(R_desc).to(torch_device)  # N, d
            R_d_desc_torch = torch.from_numpy(R_d_desc).to(torch_device)

            from .torchtools import GDMLTorchAssemble

            torch_assemble = GDMLTorchAssemble(
                J,
                tril_perms_lin,
                sig,
                R_desc_torch,
                R_d_desc_torch,
                out=K[:K_n_rows, :],
                callback=progress_callback,
            )

            # Enable data parallelism
            n_gpu = torch.cuda.device_count()
            if n_gpu > 1:
                torch_assemble = torch.nn.DataParallel(torch_assemble)
            torch_assemble.to(torch_device)

            torch_assemble.forward(torch.arange(len(J)))
            del torch_assemble

            del R_desc_torch
            del R_d_desc_torch

            stop = timeit.default_timer()

            if callback is not None:
                dur_s = stop - start
                sec_disp_str = 'took {:.1f} s'.format(dur_s) if dur_s >= 0.1 else ''
                callback(DONE, sec_disp_str=sec_disp_str)

            return K

        K = mp.RawArray('d', (K_n_rows + alloc_extra_rows) * K_n_cols)
        glob['K'], glob['K_shape'] = K, (K_n_rows + alloc_extra_rows, K_n_cols)
        glob['R_desc'], glob['R_desc_shape'] = _share_array(R_desc, 'd')
        glob['R_d_desc'], glob['R_d_desc_shape'] = _share_array(R_d_desc, 'd')

        glob['desc_func'] = desc

        start = timeit.default_timer()

        pool = None
        map_func = map
        if self._max_processes != 1 and mp.cpu_count() > 1:
            pool = Pool(
                (self._max_processes or mp.cpu_count()) - 1
            )  # exclude main process
            map_func = pool.imap_unordered

        todo, done = K_n_cols, 0
        for done_wkr in map_func(
            partial(
                _assemble_kernel_mat_wkr,
                tril_perms_lin=tril_perms_lin,
                sig=sig,
                use_E_cstr=use_E_cstr,
                exploit_sym=exploit_sym,
                cols_m_limit=cols_m_limit,
            ),
            J,
        ):
            done += done_wkr

            if callback is not None:
                callback(done, todo, newline_when_done=False)

        if pool is not None:
            pool.close()
            pool.join()  # Wait for the worker processes to terminate (to measure total runtime correctly).
            pool = None

        stop = timeit.default_timer()

        if callback is not None:
            dur_s = stop - start
            sec_disp_str = 'took {:.1f} s'.format(dur_s) if dur_s >= 0.1 else ''
            callback(DONE, sec_disp_str=sec_disp_str)

        # Release some memory.
        glob.pop('K', None)
        glob.pop('R_desc', None)
        glob.pop('R_d_desc', None)

        # return np.frombuffer(K).reshape(glob['K_shape'])
        return np.frombuffer(K).reshape((K_n_rows + alloc_extra_rows), K_n_cols)

    def draw_strat_sample(self, T, n, excl_idxs=None):
        """
        Draw sample from dataset that preserves its original distribution.

        The distribution is estimated from a histogram were the bin size is
        determined using the Freedman-Diaconis rule. This rule is designed to
        minimize the difference between the area under the empirical
        probability distribution and the area under the theoretical
        probability distribution. A reduced histogram is then constructed by
        sampling uniformly in each bin. It is intended to populate all bins
        with at least one sample in the reduced histogram, even for small
        training sizes.

        Parameters
        ----------
            T : :obj:`numpy.ndarray`
                Dataset to sample from.
            n : int
                Number of examples.
            excl_idxs : :obj:`numpy.ndarray`, optional
                Array of indices to exclude from sample.

        Returns
        -------
            :obj:`numpy.ndarray`
                Array of indices that form the sample.
        """

        if excl_idxs is None or len(excl_idxs) == 0:
            excl_idxs = None

        if n == 0:
            return np.array([], dtype=np.uint)

        if T.size == n:  # TODO: this only works if excl_idxs=None
            assert excl_idxs is None
            return np.arange(n)

        if n == 1:
            idxs_all_non_excl = np.setdiff1d(
                np.arange(T.size), excl_idxs, assume_unique=True
            )
            return np.array([np.random.choice(idxs_all_non_excl)])

        # Freedman-Diaconis rule
        h = 2 * np.subtract(*np.percentile(T, [75, 25])) / np.cbrt(n)
        n_bins = int(np.ceil((np.max(T) - np.min(T)) / h)) if h > 0 else 1
        n_bins = min(
            n_bins, int(n / 2)
        )  # Limit number of bins to half of requested subset size.

        bins = np.linspace(np.min(T), np.max(T), n_bins, endpoint=False)
        idxs = np.digitize(T, bins)

        # Exclude restricted indices.
        if excl_idxs is not None and excl_idxs.size > 0:
            idxs[excl_idxs] = n_bins + 1  # Impossible bin.

        uniq_all, cnts_all = np.unique(idxs, return_counts=True)

        # Remove restricted bin.
        if excl_idxs is not None and excl_idxs.size > 0:
            excl_bin_idx = np.where(uniq_all == n_bins + 1)
            cnts_all = np.delete(cnts_all, excl_bin_idx)
            uniq_all = np.delete(uniq_all, excl_bin_idx)

        # Compute reduced bin counts.
        reduced_cnts = np.ceil(cnts_all / np.sum(cnts_all, dtype=float) * n).astype(int)
        reduced_cnts = np.minimum(
            reduced_cnts, cnts_all
        )  # limit reduced_cnts to what is available in cnts_all

        # Reduce/increase bin counts to desired total number of points.
        reduced_cnts_delta = n - np.sum(reduced_cnts)

        while np.abs(reduced_cnts_delta) > 0:

            # How many members can we remove from an arbitrary bucket, without any bucket with more than one member going to zero?
            max_bin_reduction = np.min(reduced_cnts[np.where(reduced_cnts > 1)]) - 1

            # Generate additional bin members to fill up/drain bucket counts of subset. This array contains (repeated) bucket IDs.
            outstanding = np.random.choice(
                uniq_all,
                min(max_bin_reduction, np.abs(reduced_cnts_delta)),
                p=(reduced_cnts - 1) / np.sum(reduced_cnts - 1, dtype=float),
                replace=True,
            )
            uniq_outstanding, cnts_outstanding = np.unique(
                outstanding, return_counts=True
            )  # Aggregate bucket IDs.

            outstanding_bucket_idx = np.where(
                np.in1d(uniq_all, uniq_outstanding, assume_unique=True)
            )[
                0
            ]  # Bucket IDs to Idxs.
            reduced_cnts[outstanding_bucket_idx] += (
                np.sign(reduced_cnts_delta) * cnts_outstanding
            )
            reduced_cnts_delta = n - np.sum(reduced_cnts)

        # Draw examples for each bin.
        idxs_train = np.empty((0,), dtype=int)
        for uniq_idx, bin_cnt in zip(uniq_all, reduced_cnts):
            idx_in_bin_all = np.where(idxs.ravel() == uniq_idx)[0]
            idxs_train = np.append(
                idxs_train, np.random.choice(idx_in_bin_all, bin_cnt, replace=False)
            )

        return idxs_train
