import torch
from regression_core import SquaredExp, LazyWhite

dim = 2
kern = SquaredExp(dim=dim)
x = torch.rand(19, dim)
xx = torch.rand(37, dim)
K = kern(x=x, xx=xx)
assert torch.allclose(K, (-(x[:, None]-xx[None])**2/2)
                      .sum(dim=-1).exp())

white = LazyWhite(signal=1.0)
x = torch.rand(13, dim)
assert (white(x, x) == torch.eye(13)).all()

kern = SquaredExp(dim=dim)
white = LazyWhite(dim=dim, signal=1.0)
assert kern(x, xx, 'func').shape == white(x, xx, 'func').shape
assert kern(x, xx, 'grad').shape == white(x, xx, 'grad').shape
assert kern(x, xx, 'gradgrad').shape == white(x, xx, 'gradgrad').shape
K = white(x, operation='gradgrad')
assert (K.reshape(x.numel(), x.numel()) == (-1) *
        torch.eye(x.numel())).all()  # see NOTE 1. for -1
