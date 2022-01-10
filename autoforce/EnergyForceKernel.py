import torch
from torch.nn import Module

from util_util import iterable #, mkdir_p, safe_dirname

class EnergyForceKernel(Module):

    def __init__(self, similaritykernels):
        super().__init__()
        self.kernels = iterable(similaritykernels)
        self.name_kernels()

    def add_kernels(self, kernels):
        self.kernels = [kern for kern in self.kernels] + \
            [kern for kern in iterable(kernels)]
        self.name_kernels()

    def name_kernels(self):
        for i, kern in enumerate(self.kernels):
            kern.name = 'kern_{}'.format(i)

    @property
    def params(self):
        return [par for kern in self.kernels for par in kern.params]

    def forward(self, first, second=None, cov='energy_energy', inducing=None):
        #
        #print("\nPass here 42 in EnergyForceKernel")
        #
        sec = first if second is None else second
        #print("first = ", first)
        #print("sec = ", sec)
        if inducing is None:
            return getattr(self, cov)(first, sec)
        else:
            print("before energy_energy")
            middle = getattr(self, 'energy_energy')(inducing, inducing)
            chol, _ = jitcholesky(middle)
            invchol = chol.inverse()
            lcov, rcov = cov.split('_')
            left = getattr(self, lcov+'_energy')(first, inducing) @ invchol.t()
            if second is None and rcov == lcov:
                right = left.t()
            else:
                right = invchol @ getattr(self, 'energy_'+rcov)(inducing, sec)
            return left, right

    def energy_energy(self, first, second):
        return self.base_kerns(first, second, 'func')

    def forces_energy(self, first, second):
        return -self.base_kerns(first, second, 'leftgrad')

    def energy_forces(self, first, second):
        return -self.base_kerns(first, second, 'rightgrad')

    def forces_forces(self, first, second):
        return self.base_kerns(first, second, 'gradgrad')

    def base_kerns(self, first, second, operation):
        return torch.stack([kern(first, second, operation=operation)
                            for kern in self.kernels]).sum(dim=0)

    # diagonal elements:
    def diag(self, data, operation='energy'):
        return getattr(self, operation+'_diag')(data)

    def full_diag(self, data):
        return self.energy_forces_diag(data)

    def energy_forces_diag(self, data):
        return torch.cat([self.energy_diag(data), self.forces_diag(data)])

    def energy_diag(self, data):
        return self.base_kerns_diag(data, 'func')

    def forces_diag(self, data):
        return self.base_kerns_diag(data, 'gradgrad')

    def base_kerns_diag(self, data, operation):
        return torch.stack([kern.diag(data, operation=operation)
                            for kern in self.kernels]).sum(dim=0)

    @property
    def method_caching(self):
        return [kern.method_caching if hasattr(kern, 'method_caching') else False
                for kern in self.kernels]

    @method_caching.setter
    def method_caching(self, value):
        if hasattr(value, '__iter__'):
            val = value
        else:
            val = len(self.kernels)*[value]
        for kern, v in zip(*[self.kernels, val]):
            kern.method_caching = v

    def clear_cached(self):
        for kern in self.kernels:
            try:
                kern.cached.clear()
            except AttributeError:
                pass

    @property
    def state_args(self):
        return '[{}]'.format(', '.join([kern.state for kern in self.kernels]))

    @property
    def state(self):
        return 'EnergyForceKernel({})'.format(self.state_args)