# +
import torch
from torch.distributions import MultivariateNormal, LowRankMultivariateNormal
from collections import Counter
from scipy.optimize import minimize
from math import pi
import copy
import functools
import warnings

from regression_algebra import jitcholesky, projected_process_auxiliary_matrices_D

class ConstMean:

    def __init__(self):
        self.per_atom = torch.zeros([])
        self.unique_params = []

    def set_data(self, data):
        n = torch.as_tensor(data.natoms).view(-1)
        e = torch.stack([atoms.target_energy for atoms in data]).view(-1)
        self.per_atom = (e/n).mean()

    def __call__(self, atoms, forces=False):
        n = len(atoms)
        e = n*self.per_atom
        if forces:
            return e, torch.zeros(n, 3)
        else:
            return e


def mean_energy_per_atom_type(data):
    types, counts = zip(*[(z, c) for z, c in data.counts().items()])
    indices = {z: k for k, z in enumerate(types)}
    X = torch.zeros(len(data), len(types))
    for i, atoms in enumerate(data):
        for z, c in atoms.counts().items():
            X[i, indices[z]] = c
    N = torch.tensor(data.natoms)
    Y = data.target_energy
    mean = (Y/N).mean()
    cov = X.t()@X
    L, ridge = jitcholesky(cov)
    inv = L.cholesky_inverse()
    mu = inv@X.t()@(Y-mean*N) + mean
    weights = {z: w for z, w in zip(types, mu)}
    # max_err = (X@mu-Y).abs().max()
    return weights


class DefaultMean:

    def __init__(self, weights=None):
        self.weights = {} if weights is None else weights
        self.unique_params = []

    def set_data(self, data):
        self.weights = mean_energy_per_atom_type(data)

    def __call__(self, atoms, forces=False):
        e = torch.zeros([])
        for number, count in atoms.counts().items():
            if number in self.weights:
                e += count*self.weights[number]
        if forces:
            return e, torch.zeros(len(atoms), 3)
        else:
            return e

    def __repr__(self):
        w = {z: float(w) for z, w in self.weights.items()}
        return f"DefaultMean({w})"


class AutoMean:

    def __init__(self):
        self.weights = {}

    @property
    def unique_params(self):
        return self.weights.values()

    def sync_params(self):
        for k in sorted(self.weights.keys()):
            torch.distributed.broadcast(self.weights[k], 0)

    def set_data(self, data):
        # self._weights = mean_energy_per_atom_type(data) # unstable?
        self._weights = {z: 0. for z in data.counts().keys()}
        for z in self._weights.keys():
            if z not in self.weights:
                self.weights[z] = torch.tensor(0.)

    def __call__(self, atoms, forces=False):
        e = torch.zeros([])
        for z, count in atoms.counts().items():
            if z in self.weights:
                e += count*(self.weights[z]+self._weights[z])
        if forces:
            return e, torch.zeros(len(atoms), 3)
        else:
            return e

    def __repr__(self):
        w = {z: float(self.weights[z]+self._weights[z])
             for z in self.weights.keys()}
        return f"AutoMean({w})"


def context_setting(method):

    @functools.wraps(method)
    def wrapper(self, *args, use_caching=False, enable_grad=False, **kwargs):
        caching_status = self.gp.method_caching
        self.gp.method_caching = use_caching
        with torch.set_grad_enabled(enable_grad):
            result = method(self, *args, **kwargs)
        self.gp.method_caching = caching_status
        return result

    return wrapper



def to_0_1(x):
    return 1/x.neg().exp().add(1.)


def to_inf_inf(y):
    return (y/y.neg().add(1.)).log()


def kldiv_normal(y, sigma):
    """
    may be unstable because of often small length of y
    compared to the number of bins.
    """
    delta = sigma/10
    width = 10*sigma
    x = torch.arange(0, width, delta)
    x = torch.cat([-x.flip(0)[:-1], x])
    p = (y.view(-1)-x.view(-1, 1)).div(delta).pow(2).mul(-0.5).exp().sum(dim=1)
    nrm = torch.tensor(2*pi).sqrt()*delta
    p = (p/(y.numel()*nrm)).clamp(min=1e-8)
    q = torch.distributions.Normal(0., sigma)
    #loss = -(p*q.log_prob(x)).sum()
    loss = (p.log()-q.log_prob(x)).mul(p).sum()*delta
    return loss


def _regression(self, optimize=False, noise_f=None, max_noise=0.99, same_sigma=True, wjac=True):

    if self.ignore_forces:
        raise RuntimeError('ignore_forces is deprecated!')

    if not hasattr(self, '_noise'):
        self._noise = {}

    if noise_f is None:
        noise_f = 0.

    #
    scale = {}
    if same_sigma:
        if 'all' not in self._noise:
            self._noise['all'] = to_inf_inf(self.gp.noise.signal.detach())
        scale['all'] = self.M.diag().mean() * max_noise
    else:
        numbers = torch.tensor([x.number for x in self.X])
        zset = numbers.unique().tolist()
        for z in zset:
            if z not in self._noise:
                self._noise[z] = to_inf_inf(self.gp.noise.signal.detach())
            scale[z] = self.M.diag()[numbers == z].mean() * max_noise

    #
    L, ridge = jitcholesky(self.M)
    self.ridge = torch.as_tensor(ridge)
    self.choli = L.inverse().contiguous()
    data = self.data
    ndat = len(data)
    energies = data.target_energy
    forces = torch.cat([atoms.target_forces.view(-1) for atoms in data])
    Y = torch.cat((forces, torch.zeros(L.size(0))))

    def make_mu(with_energies=None):
        if same_sigma:
            sigma = to_0_1(self._noise['all'])*scale['all']
            sigma = sigma*torch.eye(L.shape[0])
        else:
            sigma = 0
            for z in zset:
                sigma_z = to_0_1(self._noise[z])*scale[z]
                sigma = sigma + (numbers == z).type(L.type())*sigma_z
            sigma = sigma.diag()
        if with_energies is None:
            _K = self.Kf
            _Y = Y
        else:
            _K = self.K
            _Y = torch.cat([with_energies, Y])
        #
        A = torch.cat( (_K, sigma @ L.t()) )
        Q, R = torch.linalg.qr(A)
        self.mu = (R.inverse() @ Q.t() @ _Y).contiguous()
        diff = self.Kf@self.mu - forces
        return diff

    # ------------ optimize mu ------------
    def objective_mu(x, keys, jac):
        for v, key in zip(x, keys):
            self._noise[key] = torch.tensor(v)
            if jac:
                self._noise[key].requires_grad = True
        diff = make_mu()
        loss = diff.abs().mean().sub(noise_f).pow(2)
        if jac:
            loss.backward()
            g = torch.stack([self._noise[key].grad for key in keys])
            return float(loss), g.view(-1).numpy()
        else:
            return float(loss)

    if optimize:
        # *** deprecated: extremely slow! ***
        # find approx global min on a grid;
        # if same_sigma:
        #    grid = torch.arange(.1, 5., 0.2).neg().exp()
        #    losses = {}
        #    for sig in grid:
        #        self._noise['all'] = to_inf_inf(sig)
        #        diff = make_mu()
        #        losses[sig] = diff.abs().mean().sub(noise_f).pow(2)
        #    sig = min(losses, key=losses.get)
        #    self._noise['all'] = to_inf_inf(sig)
        #    self._losses = losses
        # else:
        #    raise NotImplementedError('implement grid search!')

        # find local min
        keys = sorted(self._noise.keys())
        x0 = [float(self._noise[key]) for key in keys]
        res = minimize(objective_mu, x0=x0, jac=wjac, args=(keys, wjac))
        for v, key in zip(res.x, keys):
            self._noise[key] = torch.tensor(v)

    # make mu
    make_mu()
    if self.mu.requires_grad:
        warnings.warn('why does mu require grad?!')
        self.mu = self.mu.detach()
    self.scaled_noise = {a: float(to_0_1(b)*scale[a])
                         for a, b in self._noise.items()}

    # ------------ optimize mean ------------
    def objective_mean(w, keys, jac):
        for v, key in zip(w, keys):
            self.mean.weights[key] = torch.tensor(v)
            if jac:
                self.mean.weights[key].requires_grad = True
        mean = self.gp.mean(data, forces=False)
        diff = (mean - delta_energies)/N
        loss = diff.pow(2).mean()
        if jac:
            loss.backward()
            g = torch.stack([self.mean.weights[key].grad for key in keys])
            return float(loss), g.view(-1).numpy()
        else:
            return float(loss)

    if optimize:
        N = torch.tensor(data.natoms)
        delta_energies = energies - self.Ke@self.mu
        keys = sorted(self.mean.weights.keys())
        x0 = [float(self.mean.weights[key]) for key in keys]
        res = minimize(objective_mean, x0=x0, jac=wjac, args=(keys, wjac))
        for v, key in zip(res.x, keys):
            self.mean.weights[key] = torch.tensor(v)

    # ------------ finalize ------------
    residual = energies - self.gp.mean(data, forces=False)
    make_mu(with_energies=residual)


