clear variables; close all;

n = 2000;
p = 10;
alpha_potential = 1000;

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

fprintf('\n------- (n, p, alpha_potential) = (%d, %d, %.1f)----\n',n,p, alpha_potential);

% generate L
L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);

% intial point
X = randn(n, p);
[U, ~, V] = svd(X, 0);
X_init = U*V';
tempM2 = alpha_potential*(L\(sum(X_init.^2,2)));
tempM2 = spdiags(tempM2,0,n,n);
tempM = L + tempM2;
% XXX why need to call eigs?
[U0, ~, ~] = eigs(tempM, p, 'sm');
X0 = U0;

% Stiefel manifold
M = stiefelfactory(n,p);

% set default parameters for ARNT
opts.record = 2; % 0 for slient, 1 for outer iter. info., 2 or more for all iter. info.
opts.xtol = 1e-12;
opts.ftol = 1e-12;
opts.gtol = gtol;
opts.maxit = 500;
opts.tau = 1;
opts.usenumstab = 1;

% [~, ~, out_RGBB] = my_RGBB(X0, @PROB01_fun, M, opts, L, Lu, Ll, alpha);

% RGBB starts here
% function [x, f, out] = my_RGBB(x, fun, M, opts, varargin)
x = X0;
clear U0 X0;
fprintf("sum initial x = %18.10f\n", sum(sum(x)));

% termination rule
if ~isfield(opts, 'gtol'); opts.gtol = 1e-6;  end
if ~isfield(opts, 'xtol'); opts.xtol = 1e-6;  end
if ~isfield(opts, 'ftol'); opts.ftol = 1e-13; end

% parameters for control the linear approximation in line search,
if ~isfield(opts, 'alpha');     opts.alpha  = 1e-3;   end
if ~isfield(opts, 'rhols');     opts.rhols  = 1e-6;   end
if ~isfield(opts, 'eta');       opts.eta  = 0.2;      end
if ~isfield(opts, 'gamma');     opts.gamma  = 0.85;   end
if ~isfield(opts, 'STPEPS');    opts.STPEPS  = 1e-10; end
if ~isfield(opts, 'nt');        opts.nt  = 3;         end % 3
if ~isfield(opts, 'maxit');     opts.maxit  = 200;   end
if ~isfield(opts, 'eps');       opts.eps = 1e-14;     end
if ~isfield(opts, 'record');    opts.record = 0;      end
if ~isfield(opts, 'radius');    opts.radius = 1;      end
if isfield(opts,  'nt');         opts.nt = 5;          end

% copy parameters
gtol = opts.gtol;
xtol = opts.xtol;
ftol = opts.ftol;
maxit = opts.maxit;
rhols = opts.rhols;
eta   = opts.eta;
eps   = opts.eps;
gamma = opts.gamma;
record = opts.record;
nt = opts.nt;
alpha = opts.alpha;


% initial function value and gradient
[f, ge] = PROB01_fun(x, L, Lu, Ll, alpha_potential);
fprintf('Initial function value = %18.10f\n', f);
g = M.egrad2rgrad(x, ge);

% norm
nrmG = norm(g, 'fro');

% initial iter. information 
out.nfe = 1;
Q = 1;
Cval = f;
out.fval0 = f;

out.msg = 'exceed max iteration';
if record
    str1 = '    %6s';
    stra = ['%6s','%12s  ','%12s  ',str1,str1,str1,'   %.5s','  %.6s','\n'];
    str_head = sprintf(stra,...
        'iter', 'F','Cval', 'nrmG', 'XDiff', 'FDiff', 'nls', 'alpha');
    str1 = '  %3.2e';
    str_num = ['%4d','  %18.10e', '  %18.10e', str1,str1,str1, '  %d','  %3.2e','\n'];
end

% loop
for iter = 1:maxit
    
    xp = x;
    gp = g;
    fp = f; 
    nls = 1;
    deriv = rhols*nrmG^2; 
    d = M.lincomb(x, -1, g);
    
    % curvilinear search
    while 1
        %
        x = M.retr(xp, d, alpha);
        % Evaluate function and its gradient
        [f, ge] = PROB01_fun(x, L, Lu, Ll, alpha_potential);
        out.nfe = out.nfe + 1;
        if f <=  Cval - alpha*deriv || nls >= 5
            break
        end
        alpha = eta*alpha;
        nls = nls+1;
    end
    
    % Riemannian gradient
    g = M.egrad2rgrad(x, ge);
    
    % norm
    nrmG = M.norm(x,g);
    
    out.nrmGvec(iter) = nrmG;
    
    % difference of x
    s = x - xp;
    
    XDiff = norm(s,'inf')/alpha; % (relative Xdiff) ~ g
    FDiff = abs(f-fp)/(abs(fp)+1);
    
    % ---- record ----
    if record 
        if iter == 1
            fprintf('\n%s', str_head);
        end
        fprintf(str_num, iter, f, Cval, nrmG, XDiff, FDiff, nls, alpha);
    end
    
    % FIXME: preallocate this?
    % check stopping
    crit(iter) = FDiff;
    mcrit = mean(crit(iter-min(nt,iter)+1:iter));
    
    % ---- termination ----
    if nrmG < gtol || XDiff < xtol || FDiff < ftol
        %     if nrmG < gtol || XDiff < xtol || mcrit < ftol
        %    if nrmG < gtol
        out.msg = 'converge';
        if nrmG  < gtol, out.msg = strcat(out.msg,'_g'); end
        if XDiff < xtol, out.msg = strcat(out.msg,'_x'); end
        %         if FDiff < ftol, out.msg = strcat(out.msg,'_f'); end
        if mcrit < ftol, out.msg = strcat(out.msg,'_mf'); end
        break;
    end
    
    % difference of gradient
    if isstruct(g)
        y = matG - matGP;
    else
        y = g - gp;
    end
    
    % BB step size
    sy = abs(euclidian_iprod(s,y));
    if sy > 0
        if mod(iter,2)==0
            alpha = norm(s, 'fro')^2/sy;
        else
            alpha = sy/ norm(y, 'fro')^2;
        end
        % safeguarding on alpha
        alpha = max(min(alpha, 1e20), 1e-20);
    end
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + f)/Q;
           
end
out.XDiff = XDiff;
out.FDiff = FDiff;
out.mcrit = mcrit;
out.nrmG = nrmG;
out.fval = f;
out.iter = iter;

% Euclidean inner product
function a = euclidian_iprod(x,y)
  a = real(sum(sum(conj(x).*y)));
end
