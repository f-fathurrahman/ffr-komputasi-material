import torch
from regression_core import SquaredExp
from regression_stationary import RBF

from deprecated_regression import Covariance

x = torch.rand(23, 7)
xx = torch.rand(19, 7)
l = torch.rand(7)
s = 2.7
old = Covariance(SquaredExp(signal=s, dim=7, scale=l))
new = RBF(signal=s**2, lengthscale=l)
func = old(x, xx).allclose(new(x, xx))
leftgrad = old.leftgrad(x, xx).allclose(new.leftgrad(
    x, xx).permute(1, 0, 2).reshape(x.numel(), xx.size(0)))
rightgrad = old.rightgrad(x, xx).allclose(new.rightgrad(
    x, xx).permute(1, 2, 0).reshape(x.size(0), xx.numel()))
gradgrad = old.gradgrad(x, xx).allclose(new.gradgrad(
    x, xx).permute(2, 0, 3, 1).reshape(x.numel(), xx.numel()))
print('Squared-Exponential kernel with two different methods: \n{}\n{}\n{}\n{}'.format(
    func, leftgrad, rightgrad, gradgrad))

# test if backward works
new(x, xx).sum().backward()
