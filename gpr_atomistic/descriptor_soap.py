# +
import torch
from torch.nn import Module, Parameter
from math import factorial as fac

from descriptor_ylm import Ylm
from util_util import iterable
from regression_algebra import positive, free_form
from descriptor_func import I, Exp


class HeteroSoap(Module):

    def __init__(self, lmax, nmax, radial, numbers, atomic_unit=None, flatten=True):
        super().__init__()
        self.ylm = Ylm(lmax)
        self.nmax = nmax

        self._radial = radial
        if atomic_unit:
            self.unit = atomic_unit
        else:
            self.unit = radial.rc/3
        self.radial = Exp(-0.5*I()**2/self.unit**2)*radial
        self.numbers = sorted(iterable(numbers))
        self.species = len(self.numbers)

        one = torch.ones(lmax+1, lmax+1)
        self.Yr = 2*torch.torch.tril(one) - torch.eye(lmax+1)
        self.Yi = 2*torch.torch.triu(one, diagonal=1)

        a = torch.tensor([[1./((2*l+1)*2**(2*n+l)*fac(n)*fac(n+l))
                           for l in range(lmax+1)] for n in range(nmax+1)])
        self.nnl = (a[None]*a[:, None]).sqrt()

        self.dim = self.species**2 * (nmax+1)**2 * (lmax+1)
        self.shape = (self.species, self.species, nmax+1, nmax+1, lmax+1)
        if flatten:
            self.shape = (self.dim,)

        self.params = []
        self._state = 'atomic_unit={}, flatten={}'.format(self.unit, flatten)

        #m = torch.arange(len(self.numbers))
        #n = torch.arange(self.nmax+1)
        #i = torch.ones(m.size(0), m.size(0), n.size(0), n.size(0))
        # self.mask = (i * (m[:, None] >= m[None]).to(torch.int)[..., None, None] *
        #             (n[:, None] >= n[None]).to(torch.int)).to(torch.bool)

    @property
    def state_args(self):
        return "{}, {}, {}, {}, {}".format(self.ylm.lmax, self.nmax, self._radial.state,
                                           self.numbers, self._state)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)

    def forward(self, coo, numbers, grad=True):
        xyz = coo/self.unit
        d = xyz.pow(2).sum(dim=-1).sqrt()
        n = 2*torch.arange(self.nmax+1).type(xyz.type())
        r, dr = self.radial(self.unit*d)
        dr = self.unit*dr
        f = (r*d[None]**n[:, None])
        Y = self.ylm(xyz, grad=grad)
        if grad:
            Y, dY = Y
        ff = f[:, None, None]*Y[None]
        i = torch.arange(r.size(0))
        c = []
        for num in self.numbers:
            t = torch.index_select(ff, -1, i[numbers == num])
            c += [t.sum(dim=-1)]
        c = torch.stack(c)
        nnp = c[None, :, None, ]*c[:, None, :, None]
        p = (nnp*self.Yr).sum(dim=-1) + (nnp*self.Yi).sum(dim=-2)
        if grad:
            df = dr*d[None]**n[:, None] + r*n[:, None]*d[None]**(n[:, None]-1)
            df = df[..., None]*xyz/d[:, None]
            dc = (df[:, None, None]*Y[None, ..., None] +
                  f[:, None, None, :, None]*dY[None])
            dc = torch.stack([(numbers == num).type(r.type())[:, None] * dc
                              for num in self.numbers])
            dnnp = (c[None, :, None, ..., None, None]*dc[:, None, :, None] +
                    dc[None, :, None, ]*c[:, None, :, None, ..., None, None])
            dp = ((dnnp*self.Yr[..., None, None]).sum(dim=-3) +
                  (dnnp*self.Yi[..., None, None]).sum(dim=-4))
            p, dp = p*self.nnl, dp*self.nnl[..., None, None]/self.unit
            return p.view(*self.shape), dp.view(*self.shape, *xyz.size())
        else:
            p = p*self.nnl
            return p.view(*self.shape)


class AbsSeriesSoap(Module):

    def __init__(self, lmax, nmax, radial, unit=None):
        super().__init__()
        self.ylm = Ylm(lmax)
        self.nmax = nmax
        self.radial = radial
        if unit:
            self.unit = unit
        else:
            self.unit = radial.rc/3
        one = torch.ones(lmax+1, lmax+1)
        self.Yr = 2*torch.torch.tril(one) - torch.eye(lmax+1)
        self.Yi = 2*torch.torch.triu(one, diagonal=1)

    @property
    def state_args(self):
        return "{}, {}, {}, unit={}".format(self.ylm.lmax, self.nmax, self.radial.state, self.unit)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)

    def forward(self, coo, grad=True):
        xyz = coo/self.unit
        d = xyz.pow(2).sum(dim=-1).sqrt()
        n = 2*torch.arange(self.nmax+1).type(xyz.type())
        r, dr = self.radial(self.unit*d)
        dr = self.unit*dr
        f = (r*d[None]**n[:, None])
        Y = self.ylm(xyz, grad=grad)
        if grad:
            Y, dY = Y
        c = (f[:, None, None]*Y[None]).sum(dim=-1)
        nnp = c[None, ]*c[:, None]
        p = (nnp*self.Yr).sum(dim=-1) + (nnp*self.Yi).sum(dim=-2)
        if grad:
            df = dr*d[None]**n[:, None] + r*n[:, None]*d[None]**(n[:, None]-1)
            df = df[..., None]*xyz/d[:, None]
            dc = (df[:, None, None]*Y[None, ..., None] +
                  f[:, None, None, :, None]*dY[None])
            dnnp = (c[None, ..., None, None]*dc[:, None] +
                    dc[None, ]*c[:, None, ..., None, None])
            dp = ((dnnp*self.Yr[..., None, None]).sum(dim=-3) +
                  (dnnp*self.Yi[..., None, None]).sum(dim=-4))
            return p, dp/self.unit
        else:
            return p


class RealSeriesSoap(Module):

    def __init__(self, lmax, nmax, radial, atomic_unit=None):
        """radial: usually a cutoff function, should have an rc attr."""
        super().__init__()

        self.radial = radial
        if atomic_unit is None:
            atomic_unit = radial.rc/3
        R = Exp(-0.5*I()**2/atomic_unit**2)*radial
        self.abs = AbsSeriesSoap(lmax, nmax, R, unit=atomic_unit)

        a = torch.tensor([[1./((2*l+1)*2**(2*n+l)*fac(n)*fac(n+l))
                           for l in range(lmax+1)] for n in range(nmax+1)])
        self.nnl = (a[None]*a[:, None]).sqrt()

    def forward(self, xyz, grad=True):
        p = self.abs(xyz, grad=grad)
        if grad:
            p, q = p
            q = q*self.nnl[..., None, None]
        p = p*self.nnl

        if grad:
            return p, q
        else:
            return p

    @property
    def state_args(self):
        return "{}, {}, {}, atomic_unit={}".format(self.abs.ylm.lmax, self.abs.nmax,
                                                   self.radial.state, self.abs.unit)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)


class TailoredSoap(Module):

    def __init__(self, soap, corners=0, symm=False):
        super().__init__()
        self.soap = soap
        n = torch.arange(soap.abs.nmax+1)
        self.mask = ((n[:, None]-n[None]).abs() <=
                     soap.abs.nmax-corners)

        if not symm:
            self.mask = (self.mask & (n[:, None] >= n[None]))

        self._state_args = "corners={}, symm={}".format(corners, symm)
        self.params = []

    def forward(self, xyz, grad=True):
        p = self.soap(xyz, grad=grad)
        if grad:
            p, q = p

        p = p[self.mask].view(-1)
        if grad:
            q = q[self.mask].view(p.size(0), *xyz.size())

        if grad:
            return p, q
        else:
            return p

    @property
    def dim(self):
        return self.mask.sum()*(self.soap.abs.ylm.lmax+1)

    @property
    def state_args(self):
        return "{}, {}".format(self.soap.state, self._state_args)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)


class MultiSoap(Module):
    def __init__(self, soaps):
        super().__init__()
        self.soaps = iterable(soaps)
        self.params = [par for soap in self.soaps for par in soap.params]

    def forward(self, xyz, masks, grad=True):
        p = [soap(xyz[m], grad=grad) for soap, m in zip(*[self.soaps, masks])]
        if grad:
            p, _q = zip(*p)
            n = xyz.size(0)
            i = torch.arange(n).long()
            q = torch.cat([torch.zeros(soap.dim, n, 3).index_add(1, i[m], qq)
                           for soap, m, qq in zip(*[self.soaps, masks, _q])])
        p = torch.cat(p)

        if grad:
            return p, q
        else:
            return p

    @property
    def dim(self):
        return sum([soap.dim for soap in self.soaps])

    @property
    def state_args(self):
        return "[" + ", ".join("{}".format(soap.state) for soap in self.soaps) + "]"

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)


class ScaledSoap(Module):

    def __init__(self, soap, scales=None):
        super().__init__()
        self.soap = soap
        self.params = [par for par in soap.params]
        self.scales = scales

    @property
    def scales(self):
        return positive(self._scales)

    @scales.setter
    def scales(self, value):
        if value is None:
            v = torch.ones(self.soap.dim)
        else:
            v = torch.as_tensor(value).view(-1)
        assert (v > 0).all()
        self._scales = Parameter(free_form(v))
        self.params.append(self._scales)

    def forward(self, *args, **kwargs):
        if 'grad' in kwargs:
            grad = kwargs['grad']
        else:
            grad = True

        p = self.soap(*args, **kwargs)
        if grad:
            p, q = p
            q = q/self.scales[..., None, None]
        p = p/self.scales

        if grad:
            return p, q
        else:
            return p

    @property
    def dim(self):
        return self.soap.dim

    @property
    def state_args(self):
        return "{}, scales={}".format(self.soap.state, self.scales.data)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)


class NormalizedSoap(Module):

    def __init__(self, soap):
        super().__init__()
        self.soap = soap
        self.params = [par for par in soap.params]

    def forward(self, *args, **kwargs):
        if 'grad' in kwargs:
            grad = kwargs['grad']
        else:
            grad = True
        p = self.soap(*args, **kwargs)
        if grad:
            p, q = p

        norm = p.norm()
        if norm > 0.0:
            norm = norm + torch.finfo().eps
            p = p/norm
            if grad:
                q = q/norm
                q = q - p[..., None, None] * (p[..., None, None] * q
                                              ).sum(dim=(0))
        if grad:
            return p, q
        else:
            return p

    @property
    def dim(self):
        return self.soap.dim

    @property
    def state_args(self):
        return "{}".format(self.soap.state)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)


class SeriesSoap(Module):
    """deprecated"""

    def __init__(self, lmax, nmax, radial, unit=None, modify=None, normalize=False,
                 cutcorners=0, symm=False):
        super().__init__()
        self.abs = AbsSeriesSoap(lmax, nmax, radial, unit=unit)

        if modify:
            a = torch.tensor([[modify**(2*n+l)/((2*l+1)*2**(2*n+l)*fac(n)*fac(n+l))
                               for l in range(lmax+1)] for n in range(nmax+1)])
            self.nnl = (a[None]*a[:, None]).sqrt()
        else:
            self.nnl = torch.ones(nmax+1, nmax+1, lmax+1)

        n = torch.arange(nmax+1)
        self.mask = ((n[:, None]-n[None]).abs() <= nmax-cutcorners).byte()
        if not symm:
            self.mask = (self.mask & (n[:, None] >= n[None]).byte())

        self.normalize = normalize

        self.kwargs = 'modify={}, normalize={}, cutcorners={}, symm={}'.format(
            modify, normalize, cutcorners, symm)

        import warnings
        warnings.warn("class {} is Deprecated".format(self.__class__.__name__))

    def forward(self, xyz, grad=True):
        p = self.abs(xyz, grad=grad)
        if grad:
            p, q = p
            q = q*self.nnl[..., None, None]
        p = p*self.nnl

        p = p[self.mask].view(-1)
        if grad:
            q = q[self.mask].view(p.size(0), *xyz.size())

        if self.normalize:
            norm = p.norm()
            if norm > 0.0:
                norm = norm + torch.finfo().eps
                p = p/norm
                if grad:
                    q = q/norm
                    q = q - p[..., None, None] * (p[..., None, None] * q
                                                  ).sum(dim=(0))
        if grad:
            return p, q
        else:
            return p

    @property
    def dim(self):
        return self.mask.sum()*(self.abs.ylm.lmax+1)

    @property
    def state_args(self):
        return "{}, {}".format(self.abs.state_args, self.kwargs)

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)


def test_validity():
    import torch
    from theforce.descriptor.cutoff import PolyCut

    xyz = torch.tensor([[0.175, 0.884, -0.87, 0.354, -0.082, 3.1],
                        [-0.791, 0.116, 0.19, -0.832, 0.184, 0.],
                        [0.387, 0.761, 0.655, -0.528, 0.973, 0.]]).t()
    xyz.requires_grad = True

    target = torch.tensor([[[0.36174603, 0.39013356, 0.43448023],
                            [0.39013356, 0.42074877, 0.46857549],
                            [0.43448023, 0.46857549, 0.5218387]],

                           [[0.2906253, 0.30558356, 0.33600938],
                            [0.30558356, 0.3246583, 0.36077952],
                            [0.33600938, 0.36077952, 0.40524778]],

                           [[0.16241845, 0.18307552, 0.20443194],
                            [0.18307552, 0.22340802, 0.26811937],
                            [0.20443194, 0.26811937, 0.34109511]]])

    s = AbsSeriesSoap(2, 2, PolyCut(3.0))
    p, dp = s(xyz)
    p = p.permute(2, 0, 1)
    print('fits pre-calculated values: {}'.format(p.allclose(target)))

    p.sum().backward()
    print('fits gradients calculated by autograd: {}'.format(
        xyz.grad.allclose(dp.sum(dim=(0, 1, 2)))))

    # test with normalization turned on
    s = SeriesSoap(3, 7, PolyCut(3.0), normalize=True)
    xyz.grad *= 0
    p, dp = s(xyz)
    p.sum().backward()
    print('fits gradients calculated by autograd (normalize=True):{}'.format(
        xyz.grad.allclose(dp.sum(dim=(0)))))

    assert s.state == eval(s.state).state

    # test if works with empty tensors
    s(torch.rand(0, 3))


def test_units():
    from theforce.descriptor.cutoff import PolyCut
    xyz = torch.tensor([[0.175, 0.884, -0.87, 0.354, -0.082, 3.1],
                        [-0.791, 0.116, 0.19, -0.832, 0.184, 0.],
                        [0.387, 0.761, 0.655, -0.528, 0.973, 0.]]).t()
    xyz = xyz*3
    cutoff = 3.0*3
    xyz.requires_grad = True

    s = SeriesSoap(3, 3, PolyCut(cutoff), normalize=True)
    p, dp = s(xyz)
    p.sum().backward()
    print('grads are consistent with larger length scale: {}'.format(
        xyz.grad.allclose(dp.sum(dim=(0)))))


def test_speed(N=100):
    from theforce.descriptor.cutoff import PolyCut
    import time
    s = AbsSeriesSoap(5, 5, PolyCut(3.0))
    start = time.time()
    for _ in range(N):
        xyz = torch.rand(30, 3)
        p = s(xyz)
    finish = time.time()
    delta = (finish-start)/N
    print("speed of {}: {} sec".format(s.state, delta))


def example():
    from theforce.descriptor.cutoff import PolyCut

    lengthscale = 2.
    cutoff = 8.
    xyz = torch.tensor([[1., 0, 0], [-1., 0, 0],
                        [0, 1., 0], [0, -1., 0],
                        [0, 0, 1.], [0, 0, -1.]]) * lengthscale
    xyz.requires_grad = True
    s = SeriesSoap(2, 2, PolyCut(cutoff), normalize=True)
    p, dp = s(xyz)
    print(p)


def test_realseriessoap():
    from theforce.descriptor.cutoff import PolyCut
    xyz = torch.tensor([[0.175, 0.884, -0.87, 0.354, -0.082, 3.1],
                        [-0.791, 0.116, 0.19, -0.832, 0.184, 0.],
                        [0.387, 0.761, 0.655, -0.528, 0.973, 0.]]).t()
    xyz = xyz*3
    cutoff = 3.0*3
    xyz.requires_grad = True

    s = NormalizedSoap(TailoredSoap(RealSeriesSoap(2, 2, PolyCut(cutoff),
                                                   atomic_unit=1.5)))

    p, dp = s(xyz)
    p.sum().backward()
    test_grad = xyz.grad.allclose(dp.sum(dim=(0)))
    err_grad = (xyz.grad-dp.sum(dim=(0))).abs().max()
    print('RealSeriesSoap: grads are consistent with autograd: {} ({})'.format(
        test_grad, err_grad))
    assert eval(s.state).state == s.state


def test_multisoap():
    from theforce.descriptor.cutoff import PolyCut
    from torch import tensor
    xyz = torch.tensor([[0.175, 0.884, -0.87, 0.354, -0.082, 3.1],
                        [-0.791, 0.116, 0.19, -0.832, 0.184, 0.],
                        [0.387, 0.761, 0.655, -0.528, 0.973, 0.]]).t()
    xyz = xyz*3
    cutoff = 3.0*3
    xyz.requires_grad = True

    soaps = [TailoredSoap(RealSeriesSoap(2, 2, PolyCut(cutoff))),
             TailoredSoap(RealSeriesSoap(3, 2, PolyCut(cutoff)))]
    ms = NormalizedSoap(ScaledSoap(MultiSoap(soaps)))

    masks = [xyz[:, 0] >= 0., xyz[:, 0] < 0.]
    a, b = ms(xyz, masks)
    a.sum().backward()
    err = (xyz.grad-b.sum(dim=0)).abs().max()
    test = xyz.grad.allclose(b.sum(dim=0))
    assert ms.dim == a.size(0)
    assert ms.state == eval(ms.state).state
    print('MultiSoap: grads are consistent with autograd: {} ({})'.format(
        test, err))


def test_heterosoap():
    import torch
    from theforce.descriptor.cutoff import PolyCut

    xyz = (torch.rand(10, 3) - 0.5) * 5
    xyz.requires_grad = True
    s = HeteroSoap(7, 5, PolyCut(8.0), [10, 18], flatten=False)
    numbers = torch.tensor(4*[10]+6*[18])
    p, dp = s(xyz, numbers)
    p.sum().backward()
    print('fits gradients calculated by autograd: {}'.format(
        xyz.grad.allclose(dp.sum(dim=(0, 1, 2, 3, 4)))))

    ss = RealSeriesSoap(7, 5, PolyCut(8.0))
    pp, dpp = ss(xyz[:4])
    print('HeteroSoap == RealSeriesSoap: {}'.format(pp.allclose(p[0, 0])))
    print('HeteroSoap == RealSeriesSoap: {}'.format(
        dpp.allclose(dp[0, 0, :, :, :, :4])))
    pp, dpp = ss(xyz[4:])
    print('HeteroSoap == RealSeriesSoap: {}'.format(pp.allclose(p[1, 1])))
    print('HeteroSoap == RealSeriesSoap: {}'.format(
        dpp.allclose(dp[1, 1, :, :, :, 4:])))

    # reshape
    s = HeteroSoap(2, 2, PolyCut(8.0), [10, 18], atomic_unit=1., flatten=True)
    p, dp = s(xyz, numbers)
    print('checking dimensions: dim={}, shape={}, grad-shape={}'.format(
        s.dim, p.shape, dp.shape))


class UniversalSoap(Module):

    def __init__(self, lmax, nmax, radial, atomic_unit=None, flatten=True, normalize=True):
        super().__init__()
        self.ylm = Ylm(lmax)
        self.nmax = nmax

        self._radial = radial
        if atomic_unit:
            self.unit = atomic_unit
        else:
            self.unit = radial.rc/6
        self.radial = Exp(-0.5*I()**2/self.unit**2)*radial

        one = torch.ones(lmax+1, lmax+1)
        self.Yr = 2*torch.torch.tril(one) - torch.eye(lmax+1)
        self.Yi = 2*torch.torch.triu(one, diagonal=1)

        a = torch.tensor([[1./((2*l+1)*2**(2*n+l)*fac(n)*fac(n+l))
                           for l in range(lmax+1)] for n in range(nmax+1)])
        self.nnl = (a[None]*a[:, None]).sqrt()

        self._shape = (nmax+1, nmax+1, lmax+1)
        self.dim = (nmax+1)*(nmax+1)*(lmax+1)
        if flatten:
            self._shape = (self.dim,)
        self._size = (119, 119, *self._shape)
        self.normalize = normalize

        self.params = []
        self._state = "{}, {}, {}, atomic_unit={}, flatten={}, normalize={}".format(
            lmax, nmax, radial.state, self.unit, flatten, normalize)

    @property
    def state_args(self):
        return self._state

    @property
    def state(self):
        return self.__class__.__name__+'({})'.format(self.state_args)

    def forward(self, coo, numbers, grad=False, normalize=None, sparse_tensor=True):
        species = torch.unique(numbers, sorted=True)
        dim0 = len(species)**2
        bcasted = torch.broadcast_tensors(species[None, ], species[:, None])
        ab = torch.cat([_.reshape(1, -1) for _ in bcasted])  # alpha, beta
        xyz = coo/self.unit
        d = xyz.pow(2).sum(dim=-1).sqrt()
        n = 2*torch.arange(self.nmax+1).type(xyz.type())
        r, dr = self.radial(self.unit*d)
        dr = self.unit*dr
        f = (r*d[None]**n[:, None])
        Y = self.ylm(xyz, grad=grad)
        if grad:
            Y, dY = Y
        ff = f[:, None, None]*Y[None]
        i = torch.arange(r.size(0))
        c = []
        for num in species:
            t = torch.index_select(ff, -1, i[numbers == num])
            c += [t.sum(dim=-1)]
        c = torch.stack(c)
        nnp = c[None, :, None, ]*c[:, None, :, None]
        p = (nnp*self.Yr).sum(dim=-1) + (nnp*self.Yi).sum(dim=-2)
        if grad:
            df = dr*d[None]**n[:, None] + r*n[:, None]*d[None]**(n[:, None]-1)
            df = df[..., None]*xyz/d[:, None]
            dc = (df[:, None, None]*Y[None, ..., None] +
                  f[:, None, None, :, None]*dY[None])
            dc = torch.stack([(numbers == num).type(r.type())[:, None] * dc
                              for num in species])
            dnnp = (c[None, :, None, ..., None, None]*dc[:, None, :, None] +
                    dc[None, :, None, ]*c[:, None, :, None, ..., None, None])
            dp = ((dnnp*self.Yr[..., None, None]).sum(dim=-3) +
                  (dnnp*self.Yi[..., None, None]).sum(dim=-4))
            p, dp = p*self.nnl, dp*self.nnl[..., None, None]/self.unit
            if (normalize if normalize else self.normalize):
                norm = p.norm() + torch.finfo().eps
                p = p/norm
                dp = dp/norm
                dp = dp - p[..., None, None] * (p[..., None, None] * dp
                                                ).sum(dim=(0, 1, 2, 3, 4))
            p = p.view(dim0, *self._shape)
            dp = dp.view(dim0, *self._shape, *xyz.size())
            if sparse_tensor:
                p = torch.sparse_coo_tensor(ab, p, size=self._size)
                dp = torch.sparse_coo_tensor(
                    ab, dp, size=(*self._size, *xyz.size()))
                return p, dp
            else:
                return ab, p, self._size, dp, (*self._size, *xyz.size())
        else:
            p = p*self.nnl
            if (normalize if normalize else self.normalize):
                norm = p.norm() + torch.finfo().eps
                p = p/norm
            if sparse_tensor:
                p = torch.sparse_coo_tensor(ab, p.view(dim0, *self._shape),
                                            size=self._size)
                return p
            else:
                return ab, p.view(dim0, *self._shape), self._size


def test_UniversalSoap():
    import torch
    from theforce.descriptor.cutoff import PolyCut
    from theforce.descriptor.soap import RealSeriesSoap

    xyz = (torch.rand(10, 3) - 0.5) * 5
    xyz.requires_grad = True
    s = UniversalSoap(3, 3, PolyCut(8.0), flatten=True)
    numbers = torch.tensor(4*[10]+6*[18])
    # test grad
    p, dp = s(xyz, numbers, grad=True)
    torch.sparse.sum(p).backward()
    print('fits gradients calculated by autograd: {}'.format(
        xyz.grad.allclose(torch.sparse.sum(dp, dim=(0, 1, 2)))))

    # test non-overlapping
    numbers = torch.tensor(4*[11]+6*[19])
    pp = s(xyz, numbers, grad=False)
    print(torch.sparse.sum(p*pp).isclose(torch.tensor(0.0)))


if __name__ == '__main__' and True:
    test_validity()
    test_units()
    test_realseriessoap()
    test_multisoap()
    test_speed()
    test_heterosoap()
    test_UniversalSoap()
