clear variables; close all;

n = 2000;
p = 10;
alpha = 100;

% set tolerance for ARNT
gtol = 1e-12;

% loop      
seed = 2010;
if exist('RandStream','file')
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
else
    disp('WARNING: Not setting random seed')
    % set seed for what?
    %randrot('state', seed);
    %randn('state',seed);
end

fprintf('\n------- (n,p,alpha) = (%d, %d, %.1f)----\n',n,p, alpha);

% generate L
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);

% intial point
X = randn(n, p);
[U, ~, V] = svd(X, 0);
X_init = U*V';
tempM2 = alpha*(L\(sum(X_init.^2,2)));
tempM2 = spdiags(tempM2,0,n,n);
tempM = L + tempM2;
% XXX why need to call eigs?
[U0, ~, ~] = eigs(tempM, p, 'sm');
X0 = U0;

% Stiefel manifold
M = stiefelfactory(n,p);

% identify Eculidean gradient and Hessian
opts.grad = @PROB01_grad;
opts.fun_extra = @PROB01_fun_extra;

opts.record = 2; % 0 for slient, 1 for outer iter. info., 2 or more for all iter. info.
opts.xtol = 1e-12;
opts.ftol = 1e-12;
opts.gtol = gtol;
opts.maxit = 1; % debug
opts.tau = 1;
opts.usenumstab = 1;

t0 = tic;
[~, ~, out_RGBB] = my_RGBB(X0, @PROB01_fun, M, opts, L, Lu, Ll, alpha);
tsolve_RGBB = toc(t0);
            
% print info. in command line
fprintf('RGBB|  f: %8.6e, nrmG: %2.1e, cpu: %4.2f, OutIter: %3d, nfe: %4d,\n',...
    out_RGBB.fval, out_RGBB.nrmG, tsolve_RGBB, out_RGBB.iter, out_RGBB.nfe);
