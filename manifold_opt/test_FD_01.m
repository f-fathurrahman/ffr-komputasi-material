clear variables; close all;
xgrid = -5:0.01:5; xgrid = xgrid';
dx = xgrid(2) - xgrid(1);
N = size(xgrid, 1);
Laplacian = (1/dx^2) * gallery('tridiag', N, 1, -2, 1);

Vpot = -0.1*exp(-5*xgrid.^2);

Ham = -0.5*Laplacian + Vpot;
HamFull = full(Ham);
[eigfunc1, eigvals1] = eig(HamFull);
[eigvals1, ind] = sort(diag(eigvals1));
eigfunc1 = eigfunc1(:,ind);

opts.record = 2; % 0 for slient, 1 for outer iter. info., 2 or more for all iter. info.
opts.xtol = 1e-12;
opts.ftol = 1e-12;
opts.gtol = 1e-12;
opts.maxit = 500; % debug

Nstates = 2;
X0 = randn(N, Nstates);
[Umat, ~, Vmat] = svd(X0, 0);
X0 = Umat*Vmat';

manifold = stiefelfactory(N, Nstates);
[X, f_end, out_RGBB] = my_RGBB(X0, @PROB02_fun, manifold, opts, Ham);

