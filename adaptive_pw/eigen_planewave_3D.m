function [lambda,phi,err_post,dof,vec_k,H,M,vext,HV,u,weight] = eigen_planewave_3D(V,L, Ec, Neig, MAX_Ec)

% solve the eigenvalue problem
% (£­Delta+V_{ext})u=\lambda u
% Hu=\lamba Mu
% lambda: eigenvalues
% u: eigenfunctions
% L: width of the domain
% Ec: energy cutoff of the planewaves
% Neig: number of eigenpairs to compute
% MAX_Ec: reference grid to calculate the residual
% err_post: the ||*||^2_{H^-1} of residual
% dof: the degrees of freedom
% vec_k:the planewave vectors
% vext: the discrete potential, vext_fft is FFT of it
% HV: the matrix of potential
% weight: weight for H^{-1} norm in Fourier space
m = 4*MAX_Ec + 1; 
N = floor(Ec);
%% solve the eigenvalue problem
% initiate the degrees of freedom
% evaluate the planewave vectors
n = 0; 
for ii = -N:N
    for j = -N:N
        for k = -N:N    
            r2 = sqrt(ii^2 + j^2+k^2);
            if r2 <= Ec
                n = n+1;
            end
        end
    end
end
dof = n;
fprintf('With Ecut = %f, number of planewaves DOF = %d \n', Ec, dof);
H = zeros(dof, dof);
HV = zeros(dof, dof);
M = zeros(dof, dof);
vec_k = zeros(dof, 3);
%%
% evaluate the planewave vectors
n = 0; 
for ii = -N:N
    for j = -N:N
        for k = -N:N    
            r2 = sqrt(ii^2 + j^2+ k^2);
            if r2 <= Ec
                n=n+1;
            vec_k(n,:) = [ii, j, k];   
            end
        end
    end
end
%%
% evaluate the discrete potential
vext = zeros(m, m, m);
for ii = 1:m
    for j = 1:m
        for k=1:m
            x = -L + (ii-1)*2*L/m;
            y = -L + (j-1)*2*L/m;
            z = -L + (k-1)*2*L/m;
            vext(ii,j,k) = V(x,y,z);
        end
    end
end
%%
% FFT for vext
vext_fft = fftn(vext)/m^3;
%%
% calculate the matrix element
for ii = 1:dof
    for j = 1:dof
        if ii == j
            H(ii,j) = norm(vec_k(ii,:),2)^2 * (pi/L)^2;
            M(ii,j) = 1.0;
        end
       % external potential H_{pq}=vext_fft(p-q)
        dk = vec_k(ii,:) - vec_k(j,:);
        for a = 1:3
            if dk(a) < 0
                dk(a) = dk(a) + m;
            end
            dk(a) = dk(a) + 1;
        end
        H(ii,j) = H(ii,j) + real(vext_fft(dk(1),dk(2),dk(3)));  
        HV(ii,j) = real(vext_fft(dk(1),dk(2),dk(3))); 
    end
end
%%
% solve the eigenvalue problem Hu=\lamba Mu 
[phi, lambda] = eigs(H, M, Neig, 'SA');
lambda = diag(lambda);
%%
% compute (\Pi_n*V)*phi for future residual calculation
Vn_phi = HV * phi;
%% compute the residual
% weight for H^{-1} norm in Fourier space
weight = zeros(m, m, m);
Nx = floor(2*MAX_Ec);
for ii = -Nx:Nx
    for j = -Nx:Nx
        for k = -Nx:Nx
            r2 = sqrt(ii^2 + j^2+ k^2);
            if r2 <= 2*MAX_Ec
                if ii < 0  iik = ii + m + 1;  else  iik = ii + 1;  end
                if j < 0  jk = j + m + 1;  else  jk = j + 1;  end
                if k < 0  kk = k + m + 1;  else  kk = k + 1;  end
                weight(iik, jk, kk) = 1.0 / (1.0 + r2^2 * (pi/L)^2);
            end
        end
    end
end
%%
% evaluate the a posteriori error ||(V-V_n)\psi_n||
    err_post=zeros(Neig);
    phi_m = zeros(m, m, m, Neig);
    Vphi_m = zeros(m, m, m, Neig);
    u=zeros(m, m, m, Neig);
    for ll = 1:Neig
        % First compute (V-V_n)\psi_n in real space
        % This can be done by inverse FFT 
        N = floor(Ec);
        n = 0;
        for ii = -N:N
            for j = -N:N
                for k = -N:N
                    r2 = sqrt(ii^2 + j^2 + k^2);
                    if r2 <= Ec
                        n = n + 1;
                        if ii < 0  iik = ii + m + 1;  else  iik = ii + 1;  end
                        if j < 0  jk = j + m + 1;  else  jk = j + 1;  end
                        if k < 0  kk = k + m + 1;  else  kk = k + 1;  end                
                        phi_m(iik, jk, kk, ll) = phi(n, ll);
                        Vphi_m(iik, jk, kk, ll) = Vn_phi(n, ll);
                    end
                end
            end
        end 
        % deltaVphi = vext .* real( ifftn(phi_m) ) - real( ifftn(Vphi_m) );
        u=ifftn(phi_m(:,:,:,ll))*(m^3)/(2*L)^(3/2);
        Vu=ifftn(Vphi_m(:,:,:,ll))*(m^3)/(2*L)^(3/2);
        deltaVphi = vext .* u - Vu;
        % Transform the residual from real space back to k-space
        % to compute the H^{-1} norm
        res_fft = fftn(deltaVphi)*((2*L)^(3/2))/(m^3);
        res_weight = reshape( real(res_fft).^2 .* weight , m^3, 1);
        % Note that the error estimator is ||*||^2_{H^-1} for energy !
        err_post(ll) =err_post(ll)+norm(res_weight, 1);
    end
    return;