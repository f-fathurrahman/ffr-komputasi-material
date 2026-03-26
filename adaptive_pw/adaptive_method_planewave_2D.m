function [lambda, phi, err, Ecut, u] = adaptive_method_planewave_2D(L, Ec0, Neig, tol, MAX_Ec, MAX_iter)

% perform the adaptive algorithm to solve the eigenvalue problem by planewaves
% lambda: eigenvalues
% err:  the final a posteriori error
% L:    size of the domain
% Ec0:  initial guess of the planewave energy cutoff
% Neig: number of eigenvalues to compute
% tol:  tolerence of the accuracy
% MAX_Ec: reference grid to calculate the residual
% MAX_iter: maximum number of refinements
% Ecut: the final energy cutoff
% u: eigenvectors

% We use the relation
% (Ec_new - Ec_old) * slope = log(err) - log(tol)
% it is better to overestimate the convergence rate than underestimate it!

% initiate the computations
const = 1.0;
err_post = 1.0;
Ec = Ec0;
Ecut = zeros(MAX_iter, 1);
err = zeros(MAX_iter, 1);
iter_k = 0;
% the residual is 
% (H-\lambda_n)\psi_n = (H_n-\lambda_n)\psi_n+(V-V_n)\psi_n = (V-V_n)\psi_n
% CHECK
% V_n\psi_n should be \Pi_n V \psi_n
% compute the potential part in Hamiltonian on a reference fine grid by FFT
%%
% loop for adaptive planewave refinement
while err_post > tol
    if iter_k > MAX_iter
        fprintf('*************\n');
        fprintf('Exit since the maximun refinement steps has been reached!');
        break;
        % error('Exit since the maximun refinement steps has been reached!');
    end
    iter_k = iter_k + 1;
    fprintf('*************\n');
    fprintf('Adaptive algorithm at the %d-th step ... E_cut = %f \n', iter_k, Ec);
    % solve the eigenvalue problem
V=@V_osc_2D;
 %(2) V=@V_r4_2D;
 %(3) V=@V_cos_2D;
 %(4) V=@V_exp2_2D;
 %(5) V=@V_rneg1_2D;
 %(6) V=@V_lnr_2D;
 %(7) V=@V_gauss_2D;
    [lambda,phi,err_posto,dof,vec_k,H,M,vext,HV,u,weight] = eigen_planewave_2D(V,L, Ec, Neig, MAX_Ec);
    err_post=err_posto(Neig);
    % store the a posteriori error and the Ecut
    err(iter_k) = err_post;
    Ecut(iter_k) = Ec;
    if iter_k < 2
        % For the first step, we simply add Ec by a constant
        Ec = Ec + const;
    else
        % predict Ecut for next step
        slope = ( log(err(iter_k-1)) - log(err(iter_k)) ) / ( Ecut(iter_k) - Ecut(iter_k-1) );
        fprintf('slope = %f\n', slope);
        Ec = Ec + (log(err_post) - log(tol)) / slope ;
    end
    if err_post > tol
    fprintf('a posteriori error = %e, new Ecut for next step = %f \n', err_post, Ec);
    else 
        fprintf('a posteriori error = %e \n', err_post);
    end
end

return;