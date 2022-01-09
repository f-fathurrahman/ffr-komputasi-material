# +
import torch
from torch.nn import Parameter

from regression_kernel import Kernel
from regression_algebra import free_form, positive


class Stationary(Kernel):
    """depends only on (x-xx)/l"""

    def __init__(self, signal=1.0, lengthscale=1.0):
        super().__init__()
        self.signal = signal
        self.lengthscale = lengthscale

    @property
    def signal(self):
        return positive(self._signal)

    @signal.setter
    def signal(self, value):
        v = torch.as_tensor(value)
        assert v > 0
        self._signal = Parameter(free_form(v))
        self.params.append(self._signal)

    @property
    def lengthscale(self):
        return positive(self._lengthscale)

    @lengthscale.setter
    def lengthscale(self, value):
        v = torch.as_tensor(value).view(-1)
        assert (v > 0).all()
        self._lengthscale = Parameter(free_form(v))
        self.params.append(self._lengthscale)

    @property
    def state_args(self):
        return 'signal={}, lengthscale={}'.format(self.signal.data, self.lengthscale.data)

    def l(self, x):
        l = self.lengthscale
        while l.dim() < x.dim():
            l = l.unsqueeze(-1)
        return l

    def r(self, x, xx):
        return (x-xx)/self.l(x)

    def get_func(self, x, xx):
        return self.signal*self.get_k(self.r(x, xx))

    def get_leftgrad(self, x, xx):
        return self.signal*self.get_dk(self.r(x, xx))/self.l(x)

    def get_rightgrad(self, x, xx):
        return -self.signal*self.get_dk(self.r(x, xx))/self.l(x)

    def get_gradgrad(self, x, xx):
        l = self.l(x)
        return -self.signal*self.get_d2k(self.r(x, xx))/(l[:, None]*l[None])

    # overload the following methods
    def get_k(self, r):
        """considered to be 1 at r=0 but not a real constraint"""
        raise NotImplementedError(
            'get_k in {}'.format(self.__class__.__name__))

    def get_dk(self, r):
        raise NotImplementedError(
            'get_dk in {}'.format(self.__class__.__name__))

    def get_d2k(self, r):
        raise NotImplementedError(
            'get_d2k in {}'.format(self.__class__.__name__))


class RBF(Stationary):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_k(self, r):
        return (-0.5*(r**2).sum(dim=0)).exp()

    def get_dk(self, r):
        return -r*self.get_k(r)

    def get_d2k(self, r):
        eye = torch.eye(r.size(0))
        while eye.dim() < r.dim()+1:
            eye = eye.unsqueeze(-1)
        return (r[:, None]*r[None]-eye)*self.get_k(r)


