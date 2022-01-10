import torch
from torch.nn import Module, Parameter

class Covariance(Module):
    """
    Calculates the covariance matrix.
    A layer between stationary (base) kernels which depend
    only on r (=x-xx) and the data (x, xx).
    """

    def __init__(self, kernels):
        super().__init__()
        self.kernels = (kernels if hasattr(kernels, '__iter__')
                        else (kernels,))
        self.params = [par for kern in self.kernels for par in kern.params]

    def base_kerns(self, x=None, xx=None, operation='func'):
        return torch.stack([kern(x=x, xx=xx, operation=operation)
                            for kern in self.kernels]).sum(dim=0)

    def diag(self, x=None, operation='func'):
        if operation == 'func':
            sign = 1
        elif operation == 'gradgrad':
            sign = -1
        else:
            raise NotImplementedError('This shoud not happen!')
        return sign * torch.stack([kern.diag(x=x, operation=operation)
                                   for kern in self.kernels]).sum(dim=0)

    def forward(self, x=None, xx=None, operation='func'):
        if hasattr(self, operation):
            return getattr(self, operation)(x=x, xx=xx)
        else:
            return self.base_kerns(x=x, xx=xx, operation=operation)

    def leftgrad(self, x=None, xx=None):
        t = self.base_kerns(x=x, xx=xx, operation='grad').permute(0, 2, 1)
        return t.contiguous().view(t.size(0)*t.size(1), t.size(2))

    def rightgrad(self, x=None, xx=None):
        t = -self.base_kerns(x=x, xx=xx, operation='grad')
        return t.view(t.size(0), t.size(1)*t.size(2))

    def gradgrad(self, x=None, xx=None):
        t = -self.base_kerns(x=x, xx=xx, operation='gradgrad')
        return t.view(t.size(0)*t.size(1), t.size(2)*t.size(3))

    @property
    def state_args(self):
        return '[{}]'.format(', '.join([kern.state for kern in self.kernels]))

    @property
    def state(self):
        return 'Covariance({})'.format(self.state_args)