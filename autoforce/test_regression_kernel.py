import torch

from regression_kernel import SqD, Real, ScaledInput, Normed, DotProd, test_kernel_gradients
from regression_core import SquaredExp
from deprecated_regression import Covariance

import torch
torch.set_default_tensor_type(torch.DoubleTensor)

x = torch.rand(23, 7)
xx = torch.rand(19, 7)
old = Covariance(SquaredExp(dim=7))
new = (SqD()*Real(-0.5)).exp()
func = old(x, xx).allclose(new(x, xx))
leftgrad = old.leftgrad(x, xx).allclose(new.leftgrad(
    x, xx).permute(1, 0, 2).reshape(x.numel(), xx.size(0)))
rightgrad = old.rightgrad(x, xx).allclose(new.rightgrad(
    x, xx).permute(1, 2, 0).reshape(x.size(0), xx.numel()))
gradgrad = old.gradgrad(x, xx).allclose(new.gradgrad(
    x, xx).permute(2, 0, 3, 1).reshape(x.numel(), xx.numel()))
print('Squared-Exponential kernel with two different methods: \n{}\n{}\n{}\n{}'.format(
    func, leftgrad, rightgrad, gradgrad))

# try empty tensor
x = torch.rand(0, 7)
new(x, xx)

# test kernels gradients
kern = ScaledInput(Normed(DotProd()), scale=torch.rand(3))
print('test gradients of kernel: {}'.format(kern.state))
test_kernel_gradients(kern)
from torch import tensor
assert eval(kern.state).state == kern.state

