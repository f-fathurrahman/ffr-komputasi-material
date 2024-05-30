import numpy as np
from my_predict import GDMLPredict

def my_recov_int_const(
    model, task, R_desc=None, R_d_desc=None
):  # TODO: document e_err_inconsist return

    gdml_predict = GDMLPredict(
        model,
        max_processes=1,
        use_torch=False,
        log_level=0,
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
        print("The provided dataset contains gradients instead of force labels (flipped sign).")
        return None

    if corrcoef < 0.95:
        print("WARNING: Inconsistent energy labels detected!")
        return None

    if np.abs(e_fact - 1) > 1e-1:
        print("WARNING: Different scales in energy vs. force labels detected!")
        print(f"The integrated forces differ from the energy labels by factor ~{e_fact:.2f}")
        print("meaning that the trained model will likely fail to predict energies accurately")
        return None

    # Least squares estimate for integration constant.
    return np.sum(E_ref - E_pred) / E_ref.shape[0]

