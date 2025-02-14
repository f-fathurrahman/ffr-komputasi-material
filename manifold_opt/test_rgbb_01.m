%-------------------
function test_rgbb_01
%-------------------

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
opts.hess = @hess;
opts.grad = @grad;
opts.fun_extra = @fun_extra;

% set default parameters for ARNT
opts.record = 2; % 0 for slient, 1 for outer iter. info., 2 or more for all iter. info.
opts.xtol = 1e-12;
opts.ftol = 1e-12;
opts.gtol = gtol;
opts.maxit = 500;
opts.tau = 1;
opts.usenumstab = 1;

% run ARNT
t0 = tic;
[~, ~, out_RGBB] = RGBB(X0, @fun, M, opts);
tsolve_RGBB = toc(t0);
            
% print info. in command line
fprintf('RGBB|  f: %8.6e, nrmG: %2.1e, cpu: %4.2f, OutIter: %3d, nfe: %4d,\n',...
    out_RGBB.fval, out_RGBB.nrmG, tsolve_RGBB, out_RGBB.iter, out_RGBB.nfe);


%
% Inner functions
%
function [f,g] = fun(X,~)
    LX = L*X;
    rhoX = sum(X.^2, 2); % diag(X*X');
    tempa = Lu\(Ll\rhoX); tempa = alpha*tempa;
    
    f = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
    g = LX + bsxfun(@times,tempa,X);
end

function data = fun_extra(data)
    
    XP = data.XP;
    data.rhoX = sum(XP.^2,2);
    
end

function g = grad(X)
    rhoX = sum(X.^2, 2); % diag(X*X');
    tempa = Lu\(Ll\rhoX); tempa = alpha*tempa;
    g = L*X + bsxfun(@times,tempa,X);
end


function h = hess(X, U, data)
    
    rhoX = data.rhoX;
    rhoXdot = 2*sum(X.*U, 2);
    tempa = Lu\(Ll\rhoXdot);
    tempa = alpha*tempa;
    tempb = Lu\(Ll\rhoX);
    tempb = alpha*tempb;
    if isfield(data,'sigma')
        h = L*U + bsxfun( @times,tempa,X) + bsxfun(@times, tempb + data.sigma, U);
    else
        h = L*U + bsxfun( @times,tempa,X) + bsxfun(@times, tempb, U);
    end
end


end % function
