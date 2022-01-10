# +
import torch
from torch.nn import Module, Parameter
from regression_algebra import positive, free_form


def atleast2d(t, force2d=False):
    _t = torch.as_tensor(t)
    if _t.dim() < 2:
        return _t.view(-1, 1)
    elif _t.dim() >= 2:
        if force2d:
            return _t.view(_t.size(0), torch.tensor(_t.size()[1:]).prod())
        else:
            return _t


class Kernel(Module):
    """Kernel is a function that accepts two arguments."""

    def __init__(self):
        super().__init__()
        self.params = []

    def checkout_inputs(self, x, xx=None, diag=False):
        t = atleast2d(x, force2d=False)
        s = t.size()[1:]  # features original shape
        t = atleast2d(t, force2d=True).t().contiguous()
        if diag:
            assert xx is None
            return t, t, s
        else:
            if xx is None:
                tt = t
            else:
                tt = atleast2d(xx, force2d=True).t().contiguous()
            t, tt = torch.broadcast_tensors(t[..., None], tt[:, None])
            return t, tt, s

    # main mathods
    def forward(self, x, xx=None, diag=False, method='func'):
        t, tt, s = self.checkout_inputs(x, xx, diag)
        k = getattr(self, 'get_'+method)(t, tt)
        if k.size(-1) == 1 or (not diag and k.size(-2) == 1):
            k = torch.ones_like(t[0])*k  # TODO: use .expand instead
        return k

    def func(self, x, xx=None, diag=False):
        return self.forward(x, xx, diag, method='func')

    def leftgrad(self, x, xx=None, diag=False):
        return self.forward(x, xx, diag, method='leftgrad')

    def rightgrad(self, x, xx=None, diag=False):
        return self.forward(x, xx, diag, method='rightgrad')

    def gradgrad(self, x, xx=None, diag=False):
        return self.forward(x, xx, diag, method='gradgrad')

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)

    # Operators
    def __add__(self, other):
        return Add(self, other)

    def __sub__(self, other):
        return Sub(self, other)

    def __mul__(self, other):
        return Mul(self, other)

    def __pow__(self, other):
        return Pow(self, other)

    def pow(self, eta):
        return Pow(self, eta)

    def exp(self):
        return Exp(self)

    # overload the following methods
    @property
    def state_args(self):
        return ''

    def get_func(self, x, xx):
        """output shape: (m, n)"""
        raise NotImplementedError(
            'get_func in {}'.format(self.__class__.__name__))

    def get_leftgrad(self, x, xx):
        """output shape: (d, m, n)"""
        raise NotImplementedError(
            'get_leftgrad in {}'.format(self.__class__.__name__))

    def get_rightgrad(self, x, xx):
        """output shape: (d, m, n)"""
        raise NotImplementedError(
            'get_rightgrad in {}'.format(self.__class__.__name__))

    def get_gradgrad(self, x, xx):
        """output shape: (d, d, m, n)"""
        raise NotImplementedError(
            'get_gradgrad in {}'.format(self.__class__.__name__))


class BinaryOperator(Kernel):

    def __init__(self, a, b):
        super().__init__()
        self.a = a
        self.b = b
        self.params = a.params + b.params

    @property
    def state_args(self):
        return '{}, {}'.format(self.a.state, self.b.state)


class Add(BinaryOperator):

    def __init__(self, *args):
        super().__init__(*args)

    def get_func(self, x, xx):
        return self.a.get_func(x, xx) + self.b.get_func(x, xx)

    def get_leftgrad(self, x, xx):
        k = self.a.get_leftgrad(x, xx) + self.b.get_leftgrad(x, xx)
        return k

    def get_rightgrad(self, x, xx):
        k = self.a.get_rightgrad(x, xx) + self.b.get_rightgrad(x, xx)
        return k

    def get_gradgrad(self, x, xx):
        k = self.a.get_gradgrad(x, xx) + self.b.get_gradgrad(x, xx)
        return k


class Sub(BinaryOperator):

    def __init__(self, *args):
        super().__init__(*args)

    def get_func(self, x, xx):
        return self.a.get_func(x, xx) - self.b.get_func(x, xx)

    def get_leftgrad(self, x, xx):
        k = self.a.get_leftgrad(x, xx) - self.b.get_leftgrad(x, xx)
        return k

    def get_rightgrad(self, x, xx):
        k = self.a.get_rightgrad(x, xx) - self.b.get_rightgrad(x, xx)
        return k

    def get_gradgrad(self, x, xx):
        k = self.a.get_gradgrad(x, xx) - self.b.get_gradgrad(x, xx)
        return k


class Mul(BinaryOperator):

    def __init__(self, *args):
        super().__init__(*args)

    def get_func(self, x, xx):
        return self.a.get_func(x, xx)*self.b.get_func(x, xx)

    def get_leftgrad(self, x, xx):
        k = (self.a.get_func(x, xx)*self.b.get_leftgrad(x, xx) +
             self.b.get_func(x, xx)*self.a.get_leftgrad(x, xx))
        return k

    def get_rightgrad(self, x, xx):
        k = (self.a.get_func(x, xx)*self.b.get_rightgrad(x, xx) +
             self.b.get_func(x, xx)*self.a.get_rightgrad(x, xx))
        return k

    def get_gradgrad(self, x, xx):
        k = (self.a.get_func(x, xx)*self.b.get_gradgrad(x, xx) +
             self.b.get_func(x, xx)*self.a.get_gradgrad(x, xx) +
             self.b.get_leftgrad(x, xx)[:, None]*self.a.get_rightgrad(x, xx)[None] +
             self.a.get_leftgrad(x, xx)[:, None]*self.b.get_rightgrad(x, xx)[None])
        return k


class Pow(Kernel):

    def __init__(self, kern, eta):
        super().__init__()
        self.kern = kern
        self.eta = eta
        self.params = kern.params

    @property
    def state_args(self):
        return '{}, {}'.format(self.kern.state, self.eta)

    def get_func(self, x, xx):
        return self.kern.get_func(x, xx)**self.eta

    def get_leftgrad(self, x, xx):
        k = (self.eta*self.kern.get_func(x, xx)**(self.eta-1) *
             self.kern.get_leftgrad(x, xx))
        return k

    def get_rightgrad(self, x, xx):
        k = (self.eta*self.kern.get_func(x, xx)**(self.eta-1) *
             self.kern.get_rightgrad(x, xx))
        return k

    def get_gradgrad(self, x, xx):
        k = (self.eta*self.kern.get_func(x, xx)**(self.eta-1)*self.kern.get_gradgrad(x, xx) +
             self.eta*(self.eta-1)*self.kern.get_func(x, xx)**(self.eta-2) *
             self.kern.get_leftgrad(x, xx)[:, None]*self.kern.get_rightgrad(x, xx)[None])
        return k


class Exp(Kernel):
    def __init__(self, kern):
        super().__init__()
        self.kern = kern
        self.params = kern.params

    @property
    def state_args(self):
        return '{}'.format(self.kern.state)

    def get_func(self, x, xx):
        return self.kern.get_func(x, xx).exp()

    def get_leftgrad(self, x, xx):
        k = self.kern.get_leftgrad(x, xx)*self.kern.get_func(x, xx).exp()
        return k

    def get_rightgrad(self, x, xx):
        k = self.kern.get_rightgrad(x, xx)*self.kern.get_func(x, xx).exp()
        return k

    def get_gradgrad(self, x, xx):
        k = (self.kern.get_gradgrad(x, xx) + self.kern.get_leftgrad(x, xx)[:, None] *
             self.kern.get_rightgrad(x, xx)[None])*self.kern.get_func(x, xx).exp()
        return k


class Real(Kernel):

    def __init__(self, value):
        super().__init__()
        self.value = torch.as_tensor(value)

    @property
    def state_args(self):
        return '{}'.format(self.value)

    def get_func(self, x, xx):
        return self.value.view((x.dim()-1)*[1])

    def get_leftgrad(self, x, xx):
        return torch.zeros(x.dim()*[1])

    def get_rightgrad(self, x, xx):
        return torch.zeros(x.dim()*[1])

    def get_gradgrad(self, x, xx):
        return torch.zeros((x.dim()+1)*[1])


class Positive(Kernel):

    def __init__(self, signal=1.0, requires_grad=False):
        super().__init__()
        self.signal = signal
        self.requires_grad = requires_grad

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
    def requires_grad(self):
        return self._signal.requires_grad

    @requires_grad.setter
    def requires_grad(self, value):
        self._signal.requires_grad = value

    @property
    def state_args(self):
        return 'signal={}, requires_grad={}'.format(self.signal.data, self.requires_grad)

    def get_func(self, x, xx):
        return self.signal.view((x.dim()-1)*[1])

    def get_leftgrad(self, x, xx):
        return torch.zeros(x.dim()*[1])

    def get_rightgrad(self, x, xx):
        return torch.zeros(x.dim()*[1])

    def get_gradgrad(self, x, xx):
        return torch.zeros((x.dim()+1)*[1])


class White(Kernel):

    def __init__(self, signal=1e-3, requires_grad=False):
        super().__init__()
        self.signal = signal
        self.requires_grad = requires_grad

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
    def requires_grad(self):
        return self._signal.requires_grad

    @requires_grad.setter
    def requires_grad(self, value):
        self._signal.requires_grad = value

    @property
    def state_args(self):
        return 'signal={}, requires_grad={}'.format(self.signal.data, self.requires_grad)

    def get_func(self, x, xx):
        return self.signal*x.isclose(xx).all(dim=0).type(x.type())


class SqD(Kernel):

    def __init__(self):
        super().__init__()

    @property
    def state_args(self):
        return ''

    def get_func(self, x, xx):
        return (x-xx).pow(2).sum(dim=0)

    def get_leftgrad(self, x, xx):
        return 2*(x-xx)

    def get_rightgrad(self, x, xx):
        return -2*(x-xx)

    def get_gradgrad(self, x, xx):
        # Note: output.expand(d, d, *trail) may be needed
        d = x.size(0)
        trail = x.size()[1:]
        return -2*torch.eye(d).view(d, d, *(len(trail)*[1]))


class DotProd(Kernel):

    def __init__(self):
        super().__init__()

    @property
    def state_args(self):
        return ''

    def get_func(self, x, xx):
        return (x*xx).sum(dim=0)

    def get_leftgrad(self, x, xx):
        return torch.ones_like(x)*xx

    def get_rightgrad(self, x, xx):
        return x*torch.ones_like(xx)

    def get_gradgrad(self, x, xx):
        # Note: output.expand(d, d, *trail) may be needed
        d = x.size(0)
        trail = x.size()[1:]
        return torch.eye(d).view(d, d, *(len(trail)*[1]))


class Normed(Kernel):

    def __init__(self, kern):
        super().__init__()
        self.kern = kern
        self.params = kern.params
        self.eps = torch.finfo().eps

    @property
    def state_args(self):
        return self.kern.state

    def get_func(self, x, xx):
        n = x.norm(dim=0).clamp(min=self.eps)
        nn = xx.norm(dim=0).clamp(min=self.eps)
        return self.kern.get_func(x/n, xx/nn)

    def get_leftgrad(self, x, xx):
        n = x.norm(dim=0).clamp(min=self.eps)
        nn = xx.norm(dim=0).clamp(min=self.eps)
        y = x/n
        yy = xx/nn
        f = self.kern.get_leftgrad(y, yy)
        return f/n - (f*x).sum(dim=0)*x/n**3

    def get_rightgrad(self, x, xx):
        n = x.norm(dim=0).clamp(min=self.eps)
        nn = xx.norm(dim=0).clamp(min=self.eps)
        y = x/n
        yy = xx/nn
        f = self.kern.get_rightgrad(y, yy)
        return f/nn - (f*xx).sum(dim=0)*xx/nn**3

    def get_gradgrad(self, x, xx):
        n = x.norm(dim=0).clamp(min=self.eps)
        nn = xx.norm(dim=0).clamp(min=self.eps)
        y = x/n
        yy = xx/nn
        f = self.kern.get_gradgrad(y, yy)
        gg = (f/(n*nn) - (f*x[:, None]).sum(dim=0)*x[:, None]/(n**3*nn) -
              (f*xx[None]).sum(dim=1, keepdim=True)*xx/(n*nn**3) +
              (f*x[:, None]*xx[None]).sum(dim=(0, 1)) * x[:, None]*xx[None]/(n**3*nn**3))
        return gg


class ScaledInput(Kernel):

    def __init__(self, kern, scale=1.0):
        super().__init__()
        self.kern = kern
        self.params = kern.params
        self.scale = scale

    @property
    def state_args(self):
        return '{}, scale={}'.format(self.kern.state, self.scale.data)

    @property
    def scale(self):
        return positive(self._scale)

    @scale.setter
    def scale(self, value):
        v = torch.as_tensor(value).view(-1)
        assert (v > 0).all()
        self._scale = Parameter(free_form(v))
        self.params.append(self._scale)

    def get_func(self, x, xx):
        scale = self.scale
        while scale.dim() < x.dim():
            scale = scale[..., None]
        return self.kern.get_func(x/scale, xx/scale)

    def get_leftgrad(self, x, xx):
        scale = self.scale
        while scale.dim() < x.dim():
            scale = scale[..., None]
        return self.kern.get_leftgrad(x/scale, xx/scale)/scale

    def get_rightgrad(self, x, xx):
        scale = self.scale
        while scale.dim() < x.dim():
            scale = scale[..., None]
        return self.kern.get_rightgrad(x/scale, xx/scale)/scale

    def get_gradgrad(self, x, xx):
        scale = self.scale
        while scale.dim() < x.dim():
            scale = scale[..., None]
        return self.kern.get_gradgrad(x/scale, xx/scale)/(scale[:, None]*scale[None])


def test_kernel_gradients(kern, dim=3):

    x = torch.rand(7, dim, requires_grad=True)
    xx = torch.rand(11, dim, requires_grad=True)

    kern(x, xx).sum().backward()
    g = kern.leftgrad(x, xx).sum(dim=(2)).t()
    print(g.allclose(x.grad), (g-x.grad).abs().max().data)

    g = kern.rightgrad(x, xx).sum(dim=(1)).t()
    print(g.allclose(xx.grad), (g-xx.grad).abs().max().data)

    x.grad *= 0
    kern.rightgrad(x, xx).sum().backward()
    g = kern.gradgrad(x, xx).sum(dim=(1, 3)).t()
    print(g.allclose(x.grad), (g-x.grad).abs().max().data)

    xx.grad *= 0
    kern.leftgrad(x, xx).sum().backward()
    g = kern.gradgrad(x, xx).sum(dim=(0, 2)).t()
    print(g.allclose(xx.grad), (g-xx.grad).abs().max().data)


def example():
    polynomial = Positive(requires_grad=True) * \
        (DotProd() + Positive(1e-4, requires_grad=True))**2
    squaredexp = (SqD()*Real(-0.5)).exp()


