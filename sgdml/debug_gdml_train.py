import numpy as np

from my_sgdml.utils.desc import Desc
from my_sgdml.solvers.analytic import Analytic
from my_sgdml.utils import ui

def debug_gdml_train(
    gdml_train,
    task,
    save_progr_callback=None,
    callback=None,
):

    print(">>>> ENTER debug_gdml_train")

    task = dict(task)  # make mutable
    # Need this?

    n_train, n_atoms = task['R_train'].shape[:2]

    desc = Desc(
        n_atoms,
        max_processes=gdml_train._max_processes,
    )

    n_perms = task['perms'].shape[0]
    tril_perms = np.array([Desc.perm(p) for p in task['perms']])

    dim_i = 3 * n_atoms # not used
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

    # Descriptor is calculated here
    R = task['R_train'].reshape(n_train, -1)
    R_desc, R_d_desc = desc.from_R(
        R,
        lat_and_inv=lat_and_inv,
        callback=None
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

    max_memory_bytes = gdml_train._max_memory * 1024 ** 3

    # Memory cost of analytic solver
    est_bytes_analytic = Analytic.est_memory_requirement(n_train, n_atoms)

    # Memory overhead (solver independent)
    est_bytes_overhead = y.nbytes
    est_bytes_overhead += R.nbytes
    est_bytes_overhead += R_desc.nbytes
    est_bytes_overhead += R_d_desc.nbytes

    solver_keys = {}

    byte_to_GB = 1/(1024**3)
    print(f"Estimated required memory = {est_bytes_analytic*byte_to_GB:15.5f}")
    print(f"Estimated overhead        = {est_bytes_overhead*byte_to_GB:15.5f}")
    print(f"Max memory                = {max_memory_bytes*byte_to_GB:15.5f}")
    use_analytic_solver = (
        est_bytes_analytic + est_bytes_overhead
    ) < max_memory_bytes
    print("use_analytic_solver = ", use_analytic_solver)

    # Force using analytic solver
    if not use_analytic_solver:
        print("Forced to use analytic solver")
        use_analytic_solver = True


    gdml_train.log.info(
        'Using analytic solver (expected memory requirement: ~{})'.format(
            ui.gen_memory_str(est_bytes_analytic + est_bytes_overhead)
        )
    )

    print("Using analytic solver")
    analytic = Analytic(gdml_train, desc, callback=None)
    alphas = analytic.solve(task, R_desc, R_d_desc, tril_perms_lin, y)
    print("End of finding parameters")
    print(type(alphas))
    print("alphas.shape = ", alphas.shape)

    alphas_E = None
    alphas_F = alphas
    if task['use_E_cstr']:
        alphas_E = alphas[-n_train:]
        alphas_F = alphas[:-n_train]


    print("After finding the parameters, create model which will be returned")
    print("solver_keys = ", solver_keys)

    model = gdml_train.create_model(
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
            gdml_train._recov_int_const(model, task, R_desc=R_desc, R_d_desc=R_d_desc)
            if E_train_mean is None
            else E_train_mean
        )
        print("Recover integration constant: c = ", c)
        if c is None:
            # Something does not seem right. Turn off energy predictions for this model, only output force predictions.
            model['use_E'] = False
        else:
            model['c'] = c

    print(">>>> EXIT debug_gdml_train")

    return model


