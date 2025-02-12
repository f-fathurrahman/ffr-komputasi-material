clear variables; close all;

n = 2000; % number of basis set (or grid points)
p = 10; % number of states
alpha = 100;

L = gallery('tridiag', n, -1, 2, -1);
[Ll,Lu] = lu(L);

% intial point
X = randn(n, p);
[U, ~, V] = svd(X, 0);
X_init = U*V';
tempM2 = alpha*(L \ (sum(X_init.^2, 2)));
tempM2 = spdiags(tempM2,0,n,n);
tempM = L + tempM2;
[U0, ~, ~] = eigs(tempM, p,'sm');
X0 = U0;

X = X0; % all other variables are captured
% min 0.5*Tr(X'*L*X) + alpha/4*rho(X)'*L^{dag}*rho(X), s.t. X'*X = I_k
LX = L*X;
rhoX = sum(X.^2, 2); % diag(X*X');
% Solve Poisson equation?
tempa = Lu\(Ll\rhoX);
tempa = alpha*tempa;

% function value (total energy)
f = 0.5*sum(sum(X.*(LX))) + 1/4*(rhoX'*tempa);
% gradient (apply Hamiltonian)
g = LX + bsxfun(@times,tempa,X);
