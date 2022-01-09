from regression_algebra import *

def example_sum_packed_dim():
    sizes = torch.randint(1, 10, (100,))
    Y = [torch.ones(7, size) for size in sizes]
    Y = torch.cat(Y, dim=1)
    P = sum_packed_dim(Y, sizes.tolist())
    #print(Y.size(), P.size())
    print('sum_packed_dim works: {}'.format(P.size() == torch.Size([7, 100])))



def test_gpr_algebra(n=1000):

    # test cholesky
    K = torch.ones(n, n)
    L, ridge = jitcholesky(K)
    test_cholesky = torch.allclose(torch.mm(L, L.t()), K)
    print('Cholesky: {}'.format(test_cholesky))

    # test log_normal
    Y = torch.rand(n)
    dist = torch.distributions.MultivariateNormal(torch.zeros(n), scale_tril=L)
    test_log_normal = torch.allclose(dist.log_prob(Y), log_normal(Y, K))
    print('log_normal: {}'.format(test_log_normal))

    # test select_greedy
    X = torch.rand(100, 7)
    Z = select_greedy_simple(X, 17)
    Z = select_greedy_simple(X, 17, Z=Z)

    # test solve SVD
    A = torch.diag(torch.ones(10))
    Y = torch.linspace(0, 100, 10)
    X = solve_svd(A, Y)
    test_solve_svd = torch.allclose(X, Y)
    print('solve_svd: {}'.format(test_solve_svd))


def test_iulrf(n=100, d=7, sigma=1e-4):
    torch.set_default_tensor_type(torch.DoubleTensor)  # this is necessary!
    """keep sigma higher than 1e-4"""
    Q = torch.rand(n, 7)
    D = torch.rand(n)*sigma**2
    res = inverse_using_low_rank_factor(Q, D) @ (
        Q @ Q.t() + D.diag()) - torch.eye(n)
    print("testing inverse_using_low_rank_factor:", res.abs().max().allclose(torch.zeros(1)),
          res.abs().max())


def test_PP(n=100, d=7, sigma=1e-2):
    # full Gram matrix is X@X.t()+D.diag(), where:
    X = torch.rand(n, d)
    D = (torch.ones(n)*sigma**2)
    # sparsification:
    W = X[::10]
    M = (W@W.t())
    K = (X@W.t())
    # since D = sigma^2*I, the two methods are equivalent
    Y = torch.rand(n)
    mu1, nu1 = projected_process_auxiliary_matrices_I(K, M, Y, sigma)
    mu2, nu2, _ = projected_process_auxiliary_matrices_D(K, M, Y, D)
    a = (mu1-mu2).abs().max()
    b = (nu1-nu2).abs().max()
    print('test project process (PP): {}, {}'.format(a, b))



#example_sum_packed_dim()
#test_gpr_algebra()
#test_iulrf()
test_PP()
