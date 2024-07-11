% Setup
%
% Assumption: setup_path must be called first
%
clear all; close all; %#ok<*CLALL>

filename = 'TEMP_Si2_kpt_PBE0'; % example for non-orthogonal case
%filename = 'TEMP_SiH4_quick'; % prefix of input files

% Set up inpt defaults
S = my_inpt_defaults();

% Read .inpt file
S = read_inpt(S, filename);

% Read .ion file (not needed in this case?)
S = read_ion(S, filename);




S.temp_tol = 1e-12;

% Check the cell typ (orthogonal or non-orthogonal)
if(abs(dot(S.lat_vec(1,:),S.lat_vec(2,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(2,:),S.lat_vec(3,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(3,:),S.lat_vec(1,:))) > S.temp_tol)
    fprintf('Cell type is non-orthogonal')
    S.cell_typ = 2;
end

% Lattice vectors are stored by rows
S.lat_uvec(1,:) = S.lat_vec(1,:)/norm(S.lat_vec(1,:));
S.lat_uvec(2,:) = S.lat_vec(2,:)/norm(S.lat_vec(2,:));
S.lat_uvec(3,:) = S.lat_vec(3,:)/norm(S.lat_vec(3,:));

% Set up transformation matrices for non orthogonal cells
if S.cell_typ == 2
    % Jacobian
    S.Jacb = det(S.lat_uvec');
    assert(S.Jacb > 0.0,'Volume is negative!');

    % metric_T, Gradient and laplacian transformation matrices
    S.metric_T = S.lat_uvec * S.lat_uvec' ;
    S.metric_T(1,2) = 2*S.metric_T(1,2); 
    S.metric_T(2,3) = 2*S.metric_T(2,3); 
    S.metric_T(1,3) = 2*S.metric_T(1,3);
    S.grad_T = inv(S.lat_uvec') ;
    S.lapc_T = S.grad_T * S.grad_T' ;
    S.lapc_T(1,2) = 2*S.lapc_T(1,2); 
    S.lapc_T(2,3) = 2*S.lapc_T(2,3);
    S.lapc_T(1,3) = 2*S.lapc_T(1,3);
    %
    % Convert coordinates to Cartesian
    count_prev = 0;
    count = 0;
    for ityp = 1:S.n_typ
        if(S.IsFrac(ityp) == 0)
            S.Atm(ityp).coords = transpose(S.grad_T * transpose(S.Atm(ityp).coords));
            count = count + S.Atm(ityp).n_atm_typ;
            S.Atoms(count_prev+1:count,:) = S.Atm(ityp).coords;
            count_prev = count;
        else
            count = count + S.Atm(ityp).n_atm_typ;
            count_prev = count;
        end
    end
end





% ffr: ugh this is a bit confusing
% S.BCx is BC for one direction
% S.BC is combination of all BCx, BCy, BCz
% Value of 1 of S.BCx and S.BC means the same things
% But, for periodic S.BCx is 0.
%
% BCx = 0 -> periodic, BCx = 1 -> dirichlet
if S.BC >= 0    % if user provides BOUNDARY_CONDITION: 1-4
    if(S.BC == 1) % all Dirichlet
        S.BCx = 1; S.BCy = 1; S.BCz = 1;
    elseif(S.BC == 2)
        S.BCx = 0; S.BCy = 0; S.BCz = 0;
    elseif(S.BC == 3) % surface
        S.BCx = 0; S.BCy = 0; S.BCz = 1;
    elseif(S.BC == 4) % wire
        %S.BCx = 0; S.BCy = 1; S.BCz = 1;
        S.BCx = 1; S.BCy = 1; S.BCz = 0;
    else
        error('Boundary condition should be one among {1,2,3,4}');
    end
elseif S.BCx >= 0 && S.BCy >= 0 && S.BCz >= 0 % if user provides BCx,BCy,BCz
    n_Dirichlet = S.BCx + S.BCy + S.BCz;
    if n_Dirichlet == 0
        S.BC = 2; % Periodic BC in all 3D
    elseif n_Dirichlet == 1
        S.BC = 3; % Surface, Periodic in 2D, Dirichlet in 1D
    elseif n_Dirichlet == 2
        S.BC = 4; % Wire, Periodic in 1D, Dirichlet in 2D
    elseif n_Dirichlet == 3
        S.BC = 1; % Dirichlet in all 3D
    end
else
    % if user does not provide any BC, set default to periodic in 3D
    S.BC = 2;
    S.BCx = 0; S.BCy = 0; S.BCz = 0;
end
fprintf('Boundary conditions: [%d %d %d]\n', S.BCx, S.BCy, S.BCz)
fprintf('Boundary conditions kind: %d\n', S.BC)



% S.Nx, S.Ny, S.Nz is number of intervals now
if S.Nx > 0 && S.Ny > 0 && S.Nz > 0
    S.dx = S.L1 / S.Nx;
    S.dy = S.L2 / S.Ny;
    S.dz = S.L3 / S.Nz;
elseif S.ecut > 0
    % XXX ... This is not executed ...
    S.mesh_spacing = my_Ecut2h(S.ecut, S.FDn);
    S.Nx = max(ceil(S.L1/S.mesh_spacing), S.FDn);
    S.Ny = max(ceil(S.L2/S.mesh_spacing), S.FDn);
    S.Nz = max(ceil(S.L3/S.mesh_spacing), S.FDn);
    S.dx = S.L1 / S.Nx;
    S.dy = S.L2 / S.Ny;
    S.dz = S.L3 / S.Nz;
elseif S.mesh_spacing > 0
    S.Nx = max(ceil(S.L1/S.mesh_spacing), S.FDn);
    S.Ny = max(ceil(S.L2/S.mesh_spacing), S.FDn);
    S.Nz = max(ceil(S.L3/S.mesh_spacing), S.FDn);
    S.dx = S.L1 / S.Nx;
    S.dy = S.L2 / S.Ny;
    S.dz = S.L3 / S.Nz;
end
S.dV = S.dx * S.dy * S.dz * S.Jacb;
fprintf('Grid spacing: [%f, %f, %f]\n', S.dx, S.dy, S.dz)
fprintf('Jacobian: %f\n', S.Jacb)

% Finite-difference discretization
S.Nx = S.Nx + S.BCx;
S.Ny = S.Ny + S.BCy;
S.Nz = S.Nz + S.BCz;
S.N = S.Nx * S.Ny * S.Nz;

fprintf('Number of grid points: [%d, %d, %d]\n', S.Nx, S.Ny, S.Nz)

% FDn = 1 -> 3 point
% FDn = 2 -> 5 point
% FDn = 3 -> 7 point
% FDn = 4 -> 9 point
% FDn = 5 -> 11 point
% FDn = 6 -> 13 point


% Finite difference weights of the second derivative
FDn = S.FDn;
w2 = zeros(1,FDn+1); % why row vector?
for k = 1:FDn
    w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/(k*k*factorial(FDn-k)*factorial(FDn+k));
    w2(1) = w2(1)-2*(1/(k*k));
end
S.w2 = w2;

% Finite difference weights of the first derivative
w1 = zeros(1,FDn);
for k = 1:FDn
    w1(k+1) = ((-1)^(k+1))*(factorial(FDn)^2)/(k*factorial(FDn-k)*factorial(FDn+k));
end
S.w1 = w1;


% actual inputs:
% - S.Nx, 
% - S.dx,
% - S.w1, S.w2
% - S.FDn

% S.BC = 1; % Dirichlet in all 3D
% S.BC = 2; % Periodic BC in all 3D
% S.BC = 3; % Surface, Periodic in 2D, Dirichlet in 1D
% S.BC = 4; % Wire, Periodic in 1D, Dirichlet in 2D
%
% BCx = 0 -> periodic
% BCx = 1 -> dirichlet

% Limit to cell_typ < 3
assert(S.cell_typ < 3, "cell_type must be < 3")


Nx = S.Nx;
Ny = S.Ny;
Nz = S.Nz;
n0 = S.FDn;
w1 = S.w1;
w2 = S.w2;
dx = S.dx;
dy = S.dy;
dz = S.dz;

% D_xx laplacian in 1D
%-----------------------

% Initial number of non-zeros: including ghost nodes
nnzCount = (2 * n0 + 1) * Nx;
fprintf('Initial nnzCount = %d\n', nnzCount);

% Row and column indices and the corresponding non-zero values
% used to generate sparse matrix DL11 s.t. DL11(I(k),II(k)) = V(k)
I = zeros(nnzCount,1);
V = zeros(nnzCount,1);
II = zeros(nnzCount,1);
rowCount = 1;
count = 1;
coef_dxx = 1/dx^2;
% Find non-zero entries that use forward difference
for ii = 1:Nx
    % diagonal element
    I(count) = rowCount;
    II(count) = ii;
    V(count) = w2(1)*coef_dxx ;
    % fprintf('%d %d %d\n', count, I(count), II(count))
    count = count + 1;
    % off-diagonal elements
    for q = 1:n0
        % ii + q
        I(count) = rowCount;
        II(count) = ii + q;
        V(count) = w2(q+1)*coef_dxx;
        % fprintf('%d %d %d\n', count, I(count), II(count))
        count = count + 1;
        % ii - q
        I(count) = rowCount;
        II(count) = ii - q;
        V(count) = w2(q+1)*coef_dxx;
        % fprintf('%d %d %d\n', count, I(count), II(count))
        count = count + 1;
    end
    rowCount = rowCount + 1;
end

if S.BCx == 1 % Dirichlet BC
    % Removing outside domain entries (for periodic code this is unnecessary)
    isIn = (II >= 1) & (II <= Nx);
    S.I_11 = I(isIn);
    S.II_11 = II(isIn);
    S.V_11 = V(isIn);
elseif S.BCx == 0 % Periodic BC
    S.isOutl_11 = (II < 1);
    S.isOutr_11 = (II > Nx); % Warning: Assumed influence of only neighboring cells
    S.I_11 = I;
    S.II_11 = mod(II+(Nx-1),Nx)+1;
    S.V_11 = V;
end



% D_yy laplacian in 1D
%-----------------------

% Initial number of non-zeros: including ghost nodes
nnzCount = (2 * n0 + 1) * Ny;

% Row and column indices and the corresponding non-zero values
% used to generate sparse matrix DL22 s.t. DL22(I(k),II(k)) = V(k)
I = zeros(nnzCount,1);
V = zeros(nnzCount,1);
II = zeros(nnzCount,1);
rowCount = 1;
count = 1;
coef_dyy = 1/dy^2;

% Find non-zero entries that use forward difference
for ii = 1:Ny
    % diagonal element
    I(count) = rowCount; II(count) = ii;
    V(count) = w2(1)*coef_dyy;
    count = count + 1;
    % off-diagonal elements
    for q = 1:n0
        % ii + q
        I(count) = rowCount; II(count) = ii+q;
        V(count) = w2(q+1)*coef_dyy;
        count = count + 1;
        % ii - q
        I(count) = rowCount; II(count) = ii-q;
        V(count) = w2(q+1)*coef_dyy;
        count = count + 1;
        
    end
    rowCount = rowCount + 1;
end

if S.BCy == 1
    % Removing outside domain entries (for periodic code this is unnecessary)
    isIn = (II >= 1) & (II <= Ny);
    S.I_22 = I(isIn);
    S.II_22 = II(isIn);
    S.V_22 = V(isIn);
elseif S.BCy == 0
    S.isOutl_22 = (II < 1);
    S.isOutr_22 = (II > Ny); % Warning: Assumed influence of only neighboring cells
    S.I_22 = I;
    S.II_22 = mod(II + (Ny-1), Ny) + 1;
    S.V_22 = V;
end

% D_zz laplacian in 1D
%-----------------------

% Initial number of non-zeros: including ghost nodes
nnzCount = (2 * n0 + 1) * Nz;

% Row and column indices and the corresponding non-zero values
% used to generate sparse matrix DL33 s.t. DL33(I(k),II(k)) = V(k)
I = zeros(nnzCount,1);
V = zeros(nnzCount,1);
II = zeros(nnzCount,1);
rowCount = 1;
count = 1;
coef_dzz = 1/dz^2;

% Find non-zero entries that use forward difference
for ii = 1:Nz
    % diagonal element
    I(count) = rowCount;
    II(count) = ii;
    V(count) = w2(1)*coef_dzz ;
    count = count + 1;
    % off-diagonal elements
    for q = 1:n0
        % ii + q
        I(count) = rowCount;
        II(count) = ii+q;
        V(count) = w2(q+1)*coef_dzz;
        count = count + 1;
        % ii - q
        I(count) = rowCount;
        II(count) = ii-q;
        V(count) = w2(q+1)*coef_dzz;
        count = count + 1;
        
    end
    rowCount = rowCount + 1;
end

if S.BCz == 1
    % Removing outside domain entries (for periodic code this is unnecessary)
    isIn = (II >= 1) & (II <= Nz);
    S.I_33 = I(isIn);
    S.II_33 = II(isIn);
    S.V_33 = V(isIn);
elseif S.BCz == 0
    S.isOutl_33 = (II < 1);
    S.isOutr_33 = (II > Nz); % Warning: Assumed influence of only neighboring cells
    S.I_33 = I;
    S.II_33 = mod(II + (Nz-1), Nz) + 1;
    S.V_33 = V;
end


% Additional steps for non-orthogonal cells

if S.cell_typ == 2
    % Create 1D gradient in all directions
    %---------------------------------------
    %
    % x-direction
    %-------------
    nnz_x = 2*n0*Nx;
    G = zeros(nnz_x,1);
    R = zeros(nnz_x,1);
    A = zeros(nnz_x,1);
    rowCount = 1;
    count = 1;
    coef_dx = 1/dx;
    %
    for ii = 1:Nx
        for q = 1:n0
            % ii + q
            G(count) = rowCount;
            R(count) = ii+q;
            A(count) = w1(q+1)*coef_dx;
            count = count + 1;
            % ii - q
            G(count) = rowCount;
            R(count) = ii-q;
            A(count) = -w1(q+1)*coef_dx;
            count = count + 1;
        end
        rowCount = rowCount + 1;
    end
    %
    if S.BCx == 1
        % Removing outside domain entries (for periodic code this is unnecessary)
        isIn = (R >= 1) & (R <= Nx);
        S.I_1 = G(isIn);
        S.II_1 = R(isIn);
        S.V_1 = A(isIn);
    elseif S.BCx == 0
        S.isOutl_1 = (R < 1);
        S.isOutr_1 = (R > Nx); % Warning: Assumed influence of only neighboring cells
        S.I_1 = G;
        S.II_1 = mod(R+(Nx-1),Nx)+1;
        S.V_1 = A;
    end

    % y-direction
    %-------------
    %
    nnz_y = 2*n0*Ny;
    G = zeros(nnz_y,1);
    R = zeros(nnz_y,1);
    A = zeros(nnz_y,1);
    count =1;
    rowCount =1;
    coef_dy = 1/dy;
    %
    for jj = 1:Ny
        for q = 1:n0
            % jj + q
            G(count) = rowCount;
            R(count) = jj + q;
            A(count) = w1(q+1)*coef_dy;
            count = count + 1;
            % jj - q
            G(count) = rowCount;
            R(count) = jj - q;
            A(count) = -w1(q+1)*coef_dy;
            count = count + 1;
        end
        rowCount = rowCount + 1;
    end
    %
    if S.BCy == 1
        % Removing outside domain entries (for periodic code this is unnecessary)
        isIn = (R >= 1) & (R <= Ny);
        S.I_2 = G(isIn); S.II_2 = R(isIn); S.V_2 = A(isIn);
    elseif S.BCy == 0
        S.isOutl_2 = (R<1); S.isOutr_2 = (R>Ny);
        S.I_2 = G; S.II_2 = mod(R+(Ny-1),Ny)+1; S.V_2 = A;
    end


    % z-direction
    %-------------

    nnz_z = 2*n0*Nz;
    G = zeros(nnz_z,1);
    R = zeros(nnz_z,1);
    A = zeros(nnz_z,1);
    count =1;
    rowCount =1;
    coef_dz = 1/dz;

    for kk = 1:Nz
        for q = 1:n0
            % kk + q
            G(count) = rowCount;
            R(count) = kk+q;
            A(count) = w1(q+1)*coef_dz;
            count = count + 1;
            % kk - q
            G(count) = rowCount;
            R(count) = kk-q;
            A(count) = -w1(q+1)*coef_dz;
            count = count + 1;
        end
        rowCount = rowCount + 1;
    end

    if S.BCz == 1
        % Removing outside domain entries (for periodic code this is unnecessary)
        isIn = (R >= 1) & (R <= Nz);
        S.I_3 = G(isIn);
        S.II_3 = R(isIn);
        S.V_3 = A(isIn);
    elseif S.BCz == 0
        S.isOutl_3 = (R < 1);
        S.isOutr_3 = (R > Nz);
        S.I_3 = G;
        S.II_3 = mod(R+(Nz-1),Nz) + 1;
        S.V_3 = A;
    end

end % non-orthogonal





% Create discretized Laplacian
DL11 = sparse(S.I_11, S.II_11, S.V_11, Nx, Nx);
DG1 = sparse(S.I_1, S.II_1, S.V_1, Nx, Nx);