import torch
from torch.nn import Module

from regression_kernel import White
from EnergyForceKernel import EnergyForceKernel
from regression_gppotential import AutoMean
from util_util import iterable #, mkdir_p, safe_dirname

class GaussianProcessPotential(Module):

    def __init__(self, kernels, noise=White(signal=0.01, requires_grad=True), parametric=None):
        super().__init__()
        print("Initializing GaussianProcessPotential")
        self.kern = EnergyForceKernel(kernels)
        self.noise = noise
        self.parametric = parametric or AutoMean()

    @property
    def params(self):
        p = self.kern.params + self.noise.params
        if self.parametric is not None:
            p += self.parametric.unique_params
        return p

    @property
    def descriptors(self):
        return self.kern.kernels

    @property
    def cutoff(self):
        return max([d.cutoff for d in self.descriptors])

    def add_kernels(self, kernels):
        self.kern.add_kernels(kernels)

    @property
    def species(self):
        return [a for kern in self.kern.kernels for a in iterable(kern.a)]

    @property
    def requires_grad(self):
        return [p.requires_grad for p in self.params]

    @requires_grad.setter
    def requires_grad(self, value):
        for p in self.params:
            p.requires_grad = value

    def forward(self, data, inducing=None):
        print("Pass here 263 in GaussianProcessPotential")
        if inducing is None:
            L = torch.cat([self.kern(data, cov='energy_energy'),
                           self.kern(data, cov='forces_energy')], dim=0)
            R = torch.cat([self.kern(data, cov='energy_forces'),
                           self.kern(data, cov='forces_forces')], dim=0)
            return MultivariateNormal(torch.zeros(L.size(0)),
                                      covariance_matrix=torch.cat([L, R], dim=1) +
                                      torch.eye(L.size(0))*self.noise.signal**2)
        else:
            Q = torch.cat([self.kern(data, cov='energy_energy', inducing=inducing)[0],
                           self.kern(data, cov='forces_forces', inducing=inducing)[0]], dim=0)
            return LowRankMultivariateNormal(torch.zeros(Q.size(0)), Q, self.diagonal_ridge(data))

    def diagonal_ridge(self, data, operation='full'):
        s = self.noise.signal**2
        e_diag = torch.tensor(data.natoms, dtype=s.dtype) * s
        f_diag = torch.ones(3*sum(data.natoms)) * s
        if operation == 'energy':
            return e_diag
        elif operation == 'forces':
            return f_diag
        elif operation == 'full':
            return torch.cat([e_diag, f_diag])

    def mean(self, data, forces=True, cat=True):
        if self.parametric is None:
            if forces:
                if cat:
                    return 0
                else:
                    return 0, 0
            else:
                return 0
        else:
            e = [self.parametric(sys, forces=forces) for sys in iterable(data)]
            if forces:
                e, f = zip(*e)
                e = torch.cat([_e.view(-1) for _e in e])
                f = torch.cat(f).view(-1)
                if cat:
                    return torch.cat([e, f])
                else:
                    return e, f
            else:
                return torch.cat([_e.view(-1) for _e in e])

    def Y(self, data):
        y = torch.cat([torch.tensor([sys.target_energy for sys in data])] +
                      [sys.target_forces.view(-1) for sys in data])
        return y - self.mean(data)

    def loss(self, data, Y=None, inducing=None, logprob_loss=True, cov_loss=False):
        p = self(data, inducing=inducing)
        if hasattr(p, 'cov_factor'):
            if cov_loss:
                covariance_loss = 0.5 * ((self.kern.diag(data, 'full') - p.cov_factor.pow(2).sum(dim=-1)
                                          )/self.diagonal_ridge(data)).sum()
            else:
                covariance_loss = 0
        else:
            covariance_loss = 0
        if logprob_loss:
            lp_loss = -p.log_prob(self.Y(data) if Y is None else Y)
        else:
            lp_loss = 0
        return lp_loss + covariance_loss

    @property
    def method_caching(self):
        return self.kern.method_caching

    @method_caching.setter
    def method_caching(self, value):
        self.kern.method_caching = value

    def clear_cached(self, X=None):
        if X is None:
            self.kern.clear_cached()
        else:
            for x in iterable(X):
                if hasattr(x, 'UID'):
                    UID = x.UID()
                    for a in self.cached:
                        for b in a.values():
                            for c in list(b.keys()):
                                if UID in c:
                                    del b[c]

    @property
    def cached(self):
        return [kern.cached if hasattr(kern, 'cached') else {}
                for kern in self.kern.kernels]

    @cached.setter
    def cached(self, values):
        for kern, val in zip(*[self.kern.kernels, value]):
            kern.cached = val

    def del_cached(self):
        for kern in self.kern.kernels:
            if hasattr(kern, 'cached'):
                del kern.cached

    def attach_process_group(self, group=torch.distributed.group.WORLD):
        for kern in self.kern.kernels:
            kern.process_group = group

    def detach_process_group(self):
        for kern in self.kern.kernels:
            del kern.process_group

    @property
    def state_args(self):
        return '{}, noise={}, parametric={}'.format(self.kern.state_args, self.noise.state,
                                                    self.parametric)

    @property
    def state(self):
        return 'GaussianProcessPotential({})'.format(self.state_args)

    def __repr__(self):
        return self.state

    def to_file(self, file, flag='', mode='w'):
        from util_util import one_liner
        with open(file, mode) as f:
            f.write('\n#flag: {}\n'.format(flag))
            f.write(one_liner(self.state))
            f.write('\n')