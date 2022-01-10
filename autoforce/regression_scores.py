# +
import numpy as np
import torch


def _as_tensor(pred, target):
    p = torch.as_tensor(pred)
    t = torch.as_tensor(target)
    return p, t


def maxe(pred, target):
    p, t = _as_tensor(pred, target)
    return (p-t).abs().max()


def mae(pred, target):
    p, t = _as_tensor(pred, target)
    return (p-t).abs().mean()


def rmse(pred, target):
    p, t = _as_tensor(pred, target)
    return (p-t).pow(2).mean().sqrt()


def cd(pred, target):
    p, t = _as_tensor(pred, target)
    var1 = t.var()
    var2 = (t-p).var()
    R2 = 1-var2/var1
    return R2


coeff_of_determination = cd


def get_energy_and_forces(data):
    d_e = []
    d_f = []
    d_n = []
    for d in data:
        d_e.append(d.get_potential_energy())
        d_f.append(d.get_forces().reshape(-1))
        d_n.append(d.numbers)
    d_e = np.stack(d_e).reshape(-1)
    d_f = np.concatenate(d_f)
    d_n = np.concatenate(d_n)
    d_n = d_n.reshape(-1, 1).repeat(3, axis=1).reshape(-1)
    return d_e, d_f, d_n


