function S = my_setup_defaults(S, filename)
% After reading the input files, set up remaining default values
S.temp_tol = 1e-12;

% Density tolerance for exchange-correlation
S.xc_rhotol = 1e-14;
S.xc_magtol = 1e-8;
S.xc_sigmatol = 1e-24;

% default including gradient
S.isgradient = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exchange correlation
%   Exchange:     "nox"    none                           iexch=0
%                 "slater" Slater (alpha=2/3)             iexch=1
%                 "pbex"   Perdew-Burke-Ernzenhof exch    iexch=2
%                       options: 1 -- PBE, 2 --PBEsol, 3 -- RPBE 4 --Zhang-Yang RPBE
%                 "rPW86x"  Refitted Perdew & Wang 86     iexch=3
%                 "scanx"  SCAN exchange                  iexch=4
%   
%   Correlation:  "noc"    none                           icorr=0
%                 "pz"     Perdew-Zunger                  icorr=1 
%                 "pw"     Perdew-Wang                    icorr=2
%                 "pbec"   Perdew-Burke-Ernzenhof corr    icorr=3
%                       options: 1 -- PBE, 2 --PBEsol, 3 -- RPBE
%                 "scanc"  SCAN correlation               icorr=4
%
%   Meta-GGA:     "nom"    none                           imeta=0
%                 "scan"   SCAN-Meta-GGA                  imeta=1
%                 "rscan"  rSCAN-Meta-GGA                 imeta=1
%                 "r2scan" r2SCAN-Meta-GGA                imeta=1
%
%   van der Waals "nov"    none                           ivdw=0
%                 "vdw1"   vdW-DF1                        ivdw=1
%                 "vdw2"   vdW-DF2                        ivdw=2


% decomposition of XC, ixc = [iexch,icorr imeta ivdw]
if strcmp(S.XC, 'LDA_PW')
    S.xc = 0;
    S.ixc = [1 2 0 0];
elseif strcmp(S.XC, 'LDA_PZ')
    S.xc = 1; 
    S.ixc = [1 1 0 0];
elseif strcmp(S.XC, 'GGA_PBE')
    S.xc = 2;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
    S.isgradient = 1;
elseif strcmp(S.XC, 'GGA_PBEsol')
    S.ixc = [2 3 0 0];
    S.xc_option = [2 2];
    S.isgradient = 1;
elseif strcmp(S.XC, 'GGA_RPBE')
    S.ixc = [2 3 0 0];
    S.xc_option = [3 3];
    S.isgradient = 1;
elseif strcmp(S.XC, 'vdWDF1')
    S.xc = -102; % Zhang-Yang revPBE
    S.ixc = [2 2 0 1]; % 2+4: Zhang-Yang revPBE; 2: LDA_PW86 Correlation; 0: no kinetic energy density; 1: vdW-DF1 non-linear Correlation
    S.xc_option = [4 0];
    S.vdWDFFlag = 1;
    S.isgradient = 1;
elseif strcmp(S.XC, 'vdWDF2')
    S.xc = -108; % rPW86
    S.ixc = [3 2 0 2]; % 3: rPW86; 2: LDA_PW86 Correlation; 0: no kinetic energy density; 2: vdW-DF2 non-linear Correlation
    S.vdWDFFlag = 2;
    S.isgradient = 1;
elseif strcmp(S.XC, 'SCAN')
    S.xc = 4;
    S.ixc = [4 4 1 0]; % 4: scanx; 4: scanc; 1: need kinetic energy density; 0: no vdWDF
    S.isgradient = 1;
elseif strcmp(S.XC, 'RSCAN')
    S.xc = 4;
    S.ixc = [5 5 1 0]; % 5: rscanx; 5: rscanc; 1: need kinetic energy density; 0: no vdWDF
    S.isgradient = 1;
elseif strcmp(S.XC, 'R2SCAN')
    S.xc = 4;
    S.ixc = [6 6 1 0]; % 6: r2scanx; 6: r2scanc; 1: need kinetic energy density; 0: no vdWDF
    S.isgradient = 1;
elseif strcmp(S.XC, 'HF')
    S.xc = 40;
    S.usefock = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
    S.isgradient = 1;
elseif strcmp(S.XC, 'PBE0')
    S.xc = 41;
    S.usefock = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
    S.isgradient = 1;
elseif strcmp(S.XC, 'HSE')
    S.xc = 427;
    S.usefock = 1;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
    S.isgradient = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if S.d3Flag == 1 
    if S.xc ~= 2
        fprintf('WARNING: Cannot find D3 coefficients for this functional. DFT-D3 correction calculation canceled!\n');
        S.d3Flag  = 0;
    else
        S = set_D3_coefficients(S);
    end
end

if (S.ixc(3) == 1 && S.NLCC_flag == 1)
        error('ERROR: Currently metaGGA functionals (SCAN, R2SCAN) do not support nonlinear core correction pseudopotential.\n');
end

% calculate Nelectron
S.Nelectron = 0;
for ityp = 1:S.n_typ
    S.Nelectron = S.Nelectron + S.Atm(ityp).Z * S.Atm(ityp).n_atm_typ;
end

% Save the initial positions, since later S.Atoms will be changed
S.Atoms_init = S.Atoms;

% Cychel parameters
if (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
    % Add the folder containing all the Cychel files to the path
    addpath('../Cychel');
    count_prev = 0;
    count = 0;
    for ityp = 1:S.n_typ
        if(S.IsFrac(ityp) == 0)
            S.Atm(ityp).coords = coordinateTransformation_cychel(S,S.Atm(ityp).coords,'cart2noncart_dis');
            count = count + S.Atm(ityp).n_atm_typ;
            S.Atoms(count_prev+1:count,:) = S.Atm(ityp).coords;
            count_prev = count;
        else
            count = count + S.Atm(ityp).n_atm_typ;
            count_prev = count;
        end
    end

    % Find minimum and maximum radial coordinate among all atoms
    S.xmin_at = min(S.Atoms(:,1));
    S.xmax_at = max(S.Atoms(:,1));

    % Find vacuum
    S.xvac = (S.L1 - (S.xmax_at - S.xmin_at))/2;

    % Find inner and outer radius
    S.xin  = S.xmin_at - S.xvac;
    S.xout = S.xmax_at + S.xvac;
    fprintf('Cychel radial direction parameters: \n');
    fprintf('Vacuum %f, Inner radius %f, Outer radius %f\n',S.xvac,S.xin,S.xout);

    % Rotational matrix
    theta1 = S.L2;
    S.RotM1 = [cos(theta1),-sin(theta1),0; sin(theta1),cos(theta1),0; 0 0 1];
    theta2 = S.alph*S.L3;
    S.RotM2 = [cos(theta2),-sin(theta2),0; sin(theta2),cos(theta2),0; 0 0 1];
end

% Check the cell typ (orthogonal or non-orthogonal)
if(abs(dot(S.lat_vec(1,:),S.lat_vec(2,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(2,:),S.lat_vec(3,:))) > S.temp_tol ||...
   abs(dot(S.lat_vec(3,:),S.lat_vec(1,:))) > S.temp_tol)
    S.cell_typ = 2;
end

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

% Brillouin-Zone Sampling
S = my_generate_kpts(S);

% check spin-orbit coupling
for ityp = 1:S.n_typ
    if S.Atm(ityp).pspsoc == 1
        S.SOC_flag = 1;        
        break;
    end
end

% no-spin polarized calculation
if S.spin_typ == 0
    S.nspin = 1;
    S.nspden = 1;
    if S.SOC_flag == 1
        S.nspinor = 2;
    end
% collinear polarized calculation
elseif S.spin_typ == 1
    if S.SOC_flag == 1
        error('ERROR: Collinear spin could not be used with SOC, please use non-collinear spin (SPIN_TYP: 2)!');
    end
    S.nspin = 2;
    S.nspinor = 2;
    S.nspden = 2;
% non-collinear polarized calculation
elseif S.spin_typ == 2
    S.nspin = 1;
    S.nspinor = 2;
    S.nspden = 4;
end
fprintf(' nspin = %d, nspinor = %d, nspden = %d\n', S.nspin, S.nspinor, S.nspden);

S.occfac = 2/S.nspinor;
S.nspinor_eig = S.nspinor/S.nspin;

% Provide default spin if not provided
if S.spin_typ == 1
    rng('default');
    for ityp = 1:S.n_typ
        if(S.IsSpin(ityp) == 0)
            S.Atm(ityp).mag(:,3) = -S.Atm(ityp).Z + 2 * S.Atm(ityp).Z * rand(S.Atm(ityp).n_atm_typ,1);
        end
    end
elseif S.spin_typ == 2
    rng('default');
    for ityp = 1:S.n_typ
        if(S.IsSpin(ityp) == 0)
            S.Atm(ityp).mag = -S.Atm(ityp).Z + 2 * S.Atm(ityp).Z * rand(S.Atm(ityp).n_atm_typ,3);
        end
    end
end

% check magnetization
if S.spin_typ == 1
    for ityp = 1:S.n_typ
        for i = 1:S.Atm(ityp).n_atm_typ
            if S.Atm(ityp).mag(i,1) || S.Atm(ityp).mag(i,2)
                error('ERROR: For collinear spin, the initial spin on x and y direction should be 0.')
            end
            if abs(S.Atm(ityp).mag(i,3)) > S.Atm(ityp).Z
                fprintf("WARNING: For atom type order %d, index %d, the initial magnetization is larger than zion %d.\n", ityp,i,S.Atm(ityp).Z);
            end
        end
    end
elseif S.spin_typ == 2
    for ityp = 1:S.n_typ
        for i = 1:S.Atm(ityp).n_atm_typ
            mag = sqrt(sum(S.Atm(ityp).mag(i,:).*S.Atm(ityp).mag(i,:),2));
            if mag > S.Atm(ityp).Z
                fprintf("WARNING: For atom type order %d, index %d, the initial magnetization is larger than zion %d.\n", ityp,i,S.Atm(ityp).Z);
            end
        end
    end
end

% set up default smearing if not provided
if S.bet < 0
    if S.elec_T_type == 0 % fermi-dirac 
        % The electronic temperature corresponding to 0.1 eV is 1160.452211 K
        S.bet = 27.21138602 / 0.1; % smearing = 0.1 eV = 0.00367493225 Ha, Beta := 1 / smearing
    elseif S.elec_T_type == 1 % gaussian smearing
        % The electronic temperature corresponding to 0.2 eV is 2320.904422 K
        S.bet = 27.21138602 / 0.2; % smearing = 0.2 eV = 0.00734986450 Ha, Beta := 1 / smearing
    end
    S.Temp = 1./(3.166810501187400e-06 * S.bet); 
end

% BCx = 0 -> periodic, BCx = 1 -> dirichlet
if S.BC >= 0    % if user provides BOUNDARY_CONDITION: 1-4
    if(S.BC == 1)
        S.BCx = 1; S.BCy = 1; S.BCz = 1;
    elseif(S.BC == 2)
        S.BCx = 0; S.BCy = 0; S.BCz = 0;
    elseif(S.BC == 3)
        S.BCx = 0; S.BCy = 0; S.BCz = 1;
    elseif(S.BC == 4)
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

L1 = S.L1; L2 = S.L2; L3 = S.L3;

% S.Nx, S.Ny, S.Nz is number of intervals now
if S.Nx > 0 && S.Ny > 0 && S.Nz > 0
    S.dx = S.L1 / S.Nx;
    S.dy = S.L2 / S.Ny;
    S.dz = S.L3 / S.Nz;
elseif S.ecut > 0
    S.mesh_spacing = Ecut2h(S.ecut, S.FDn);
    S.Nx = max(ceil(S.L1/S.mesh_spacing),S.FDn);
    S.Ny = max(ceil(S.L2/S.mesh_spacing),S.FDn);
    S.Nz = max(ceil(S.L3/S.mesh_spacing),S.FDn);
    S.dx = S.L1 / S.Nx;
    S.dy = S.L2 / S.Ny;
    S.dz = S.L3 / S.Nz;
elseif S.mesh_spacing > 0
    S.Nx = max(ceil(S.L1/S.mesh_spacing),S.FDn);
    S.Ny = max(ceil(S.L2/S.mesh_spacing),S.FDn);
    S.Nz = max(ceil(S.L3/S.mesh_spacing),S.FDn);
    S.dx = S.L1 / S.Nx;
    S.dy = S.L2 / S.Ny;
    S.dz = S.L3 / S.Nz;
end

dx = S.dx; dy = S.dy; dz = S.dz;
S.dV = S.dx * S.dy * S.dz * S.Jacb;

% Finite-difference discretization
S.Nx = S.Nx + S.BCx;
S.Ny = S.Ny + S.BCy;
S.Nz = S.Nz + S.BCz;
Nx = S.Nx; Ny = S.Ny; Nz = S.Nz;
S.N = S.Nx * S.Ny * S.Nz;

% map atom positions back to domain for periodic domain    
% in x direction    
isAtomOutx = sum(S.Atoms(:,1) < 0 | S.Atoms(:,1) >= S.L1) > 0;    
if (isAtomOutx)    
    if S.BCx == 0    
        S.Atoms(:,1) = mod(S.Atoms(:,1), S.L1);    
        for ityp = 1:S.n_typ    
            S.Atm(ityp).coords(:,1) = mod(S.Atm(ityp).coords(:,1),S.L1);    
        end    
        fprintf(' WARNING: mapped atom position back to domain in x dir!\n');    
    else    
        error('Atom out of domain in x dir!');    
    end    
end    
% in y direction    
isAtomOuty = sum(S.Atoms(:,2) < 0 | S.Atoms(:,2) >= S.L2) > 0;    
if (isAtomOuty)    
    if S.BCy == 0    
        S.Atoms(:,2) = mod(S.Atoms(:,2), S.L2);    
        for ityp = 1:S.n_typ    
            S.Atm(ityp).coords(:,2) = mod(S.Atm(ityp).coords(:,2),S.L2);    
        end    
        fprintf(' WARNING: mapped atom position back to domain in y dir!\n');    
    else    
        error('Atom out of domain in y dir!');    
    end    
end    
% in z direction    
isAtomOutz = sum(S.Atoms(:,3) < 0 | S.Atoms(:,3) >= S.L3) > 0;    
if (isAtomOutz)    
    if S.BCz == 0    
        S.Atoms(:,3) = mod(S.Atoms(:,3), S.L3);    
        for ityp = 1:S.n_typ    
            S.Atm(ityp).coords(:,3) = mod(S.Atm(ityp).coords(:,3),S.L3);    
        end    
        fprintf(' WARNING: mapped atom position back to domain in z dir!\n');    
    else    
        error('Atom out of domain in z dir!');    
    end    
end    

if (isAtomOutx || isAtomOuty || isAtomOutz)    
    fprintf(' After mapping\n COORD:\n');    
    disp(S.Atoms);    
end

% Finite difference weights of the second derivative
FDn = S.FDn;
w2 = zeros(1,FDn+1) ;
for k=1:FDn
    w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
        (k*k*factorial(FDn-k)*factorial(FDn+k));
    w2(1) = w2(1)-2*(1/(k*k));
end
S.w2 = w2;

% Finite difference weights of the first derivative
w1 = zeros(1,FDn) ;
for k=1:FDn
    w1(k+1) = ((-1)^(k+1))*(factorial(FDn)^2)/...
        (k*factorial(FDn-k)*factorial(FDn+k));
end
S.w1 = w1;

% Weights for spatial integration over domain
% S.W = IntgWts(S.Nx,S.Ny,S.Nz,S.BCx,S.BCy,S.BCz,S.xin,S);
if S.cell_typ == 1 || S.cell_typ == 2
    S.W = ones(S.N,1) * (S.dx*S.dy*S.dz*S.Jacb);
else
    S.W = IntgWts(S.Nx,S.Ny,S.Nz,S.BCx,S.BCy,S.BCz,S.xin,S);
end

% assert(abs(sum(S.W)-pi*(S.xout*S.xout-S.xin*S.xin)*S.L3/25)<1e-6, ...
%      'Incorrect weights for spatial integration!');

% Create spherical harmonics for poisson solve for isolated clusters
if (S.BCx == 1 && S.BCy == 1 && S.BCz == 1)
    % Calculate Spherical Harmonics with origin shifted to the center of the domain
    xx_aug = (0-FDn:Nx+FDn-1)*dx;% - L1/2;
    yy_aug = (0-FDn:Ny+FDn-1)*dy;% - L2/2;
    zz_aug = (0-FDn:Nz+FDn-1)*dz;% - L3/2;
    [XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D] = ndgrid(xx_aug,yy_aug,zz_aug);
    % Find distances
    RR_AUG_3D = calculateDistance(XX_AUG_3D,YY_AUG_3D,ZZ_AUG_3D,L1/2,L2/2,L3/2,S);

    XX_AUG = reshape(XX_AUG_3D,[],1);
    YY_AUG = reshape(YY_AUG_3D,[],1);
    ZZ_AUG = reshape(ZZ_AUG_3D,[],1);
    RR_AUG = reshape(RR_AUG_3D,[],1);

    S.RR_AUG = RR_AUG;
    S.RR_AUG_3D = RR_AUG_3D;

    in_flag = ones(Nx+2*FDn,Ny+2*FDn,Nz+2*FDn);
    in_flag(1:FDn,:,:) = 0;
    in_flag(:,1:FDn,:) = 0;
    in_flag(:,:,1:FDn) = 0;
    in_flag(FDn+Nx+1:end,:,:) = 0;
    in_flag(:,FDn+Ny+1:end,:) = 0;
    in_flag(:,:,FDn+Nz+1:end) = 0;
    isIn = (in_flag~=0);

    S.isIn = isIn;
    S.RR = RR_AUG(isIn);

    pos_node_cart = coordinateTransformation(S,[XX_AUG,YY_AUG,ZZ_AUG],'noncart2cart_dis');
    pos_atm_cart = coordinateTransformation(S,[L1/2,L2/2,L3/2],'noncart2cart_dis');
    XX_AUG = bsxfun(@minus,pos_node_cart(:,1),pos_atm_cart(:,1));
    YY_AUG = bsxfun(@minus,pos_node_cart(:,2),pos_atm_cart(:,2));
    ZZ_AUG = bsxfun(@minus,pos_node_cart(:,3),pos_atm_cart(:,3));

    l_cut = 6;
    SH = repmat(struct([]),l_cut+1,1);
    for l = 0:l_cut
        for m = -l:l
            SH(l+1).Ylm_AUG(:,m+l+1) = sphericalHarmonics(XX_AUG,YY_AUG,ZZ_AUG,l,m,'real');
            Ylm_AUG_TEMP = SH(l+1).Ylm_AUG(:,m+l+1);
            SH(l+1).Ylm(:,m+l+1) = Ylm_AUG_TEMP(isIn);
            % SH(l+1).Ylm(:,m+l+1) = sphericalHarmonics(XX,YY,ZZ,l,m,'real');
        end
    end

    S.l_cut = l_cut;
    S.SH = SH;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default values to the same default as SPARC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first find effective mesh size
if S.cell_typ < 3
    dx2_inv = 1/(S.dx * S.dx);
    dy2_inv = 1/(S.dy * S.dy);
    dz2_inv = 1/(S.dz * S.dz);
    h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
    dx2_inv = 1/(S.dx * S.dx);
    dy2_inv = 1/(((S.xin+S.xout)/2)*S.dy)^2;
    dz2_inv = 1/(S.dz * S.dz);
    h_eff = sqrt(3.0 / (dx2_inv + dy2_inv + dz2_inv));
end
% find npl
if S.npl < 0
    fprintf('## Chebyshev polynomial degree not provided, finding npl ...\n');
    S.npl = my_mesh2chebdegree(h_eff);
    fprintf('## Based on the mesh size, npl is set to: %d\n',S.npl);
end

% Nev
if S.Nev < 0
    fprintf('## Number of states not provided, finding Nev ...\n');
    S.Nev = S.nspinor_eig*(floor(S.Nelectron / 2) * 1.2 + 5); 
    S.Nev = round(S.Nev);
    fprintf('## Based on the number of electrons, Nev is set to: %d\n',S.Nev);
end

% SCF_tol
if S.SCF_tol < 0    
    if S.MDFlag 
        % in case of MD, using 1E-3 Ha/Bohr force accuracy as target
        target_force_accuracy = 1E-3;
        a = 1.025;
        b = 1.368;
        S.SCF_tol = exp((log(target_force_accuracy) - b)/a);
    elseif S.RelaxFlag
        % in case of relaxation, using TOL_RELAX/5 force accuracy as target
        target_force_accuracy = S.TOL_RELAX/5;
        a = 1.025;
        b = 1.468;
        S.SCF_tol = exp((log(target_force_accuracy) - b)/a);
    else
        % in case of single point calculation, using 1E-5 Ha/atom energy accuracy as target
        target_force_accuracy = -1.0;
        target_energy_accuracy = -1.0;
        if S.accuracy_level >= 0
            target_force_accuracy = 10^(S.accuracy_level + 1);
        elseif S.target_force_accuracy > 0
            target_force_accuracy = S.target_force_accuracy;
        elseif S.target_energy_accuracy > 0
            target_energy_accuracy = S.target_energy_accuracy;
        end
        
        % if none of the accuracy levels are specified, set energy_accuracy to
        % 1e-5
        if target_force_accuracy < 0  && target_energy_accuracy < 0 
            target_energy_accuracy = 1e-5;
        end
        
        % choose SCF tol based on the desired accuracy
        if target_energy_accuracy > 0
            a = 1.502;
            b = 1.165;
            S.SCF_tol = exp((log(target_energy_accuracy) - b)/a);
        elseif target_force_accuracy > 0
            a = 1.025;
            b = 1.368;
            S.SCF_tol = exp((log(target_force_accuracy) - b)/a);
        end
    end
    fprintf('## Based on the desired accuracy, SCF_tol is set to: %.3e\n',S.SCF_tol);
end

% poisson_tol
if S.poisson_tol < 0
    fprintf('## Poisson tolerance not provided, choosing poisson_tol ...\n')
    S.poisson_tol = S.SCF_tol * 0.01; 
    fprintf('## poisson_tol is set to: %.3e\n',S.poisson_tol);
end

% pseudocharge_tol
if S.pseudocharge_tol < 0
    fprintf('## Pseudocharge tolerance not provided, choosing pseudocharge_tol ...\n')
    S.pseudocharge_tol = S.SCF_tol * 0.001;
    fprintf('## pseudocharge_tol is set to: %.3e\n',S.pseudocharge_tol);
end

% default Kerker tolerance
if S.precond_tol < 0
    S.precond_tol = h_eff * h_eff * 0.001;
end

% mixing parameter for simple mixing
if S.MixingParameterSimple < 0
    S.MixingParameterSimple = S.MixingParameter;
end

% set default mixing parameter for magnetization density to the same as mixing
% parameter for total density/potential
if S.MixingParameterMag < 0.0
    S.MixingParameterMag = S.MixingParameter;
end

% set default simple (linear) mixing parameter for magnetization density to be the
% same as for pulay mixing
if S.MixingParameterSimpleMag < 0.0
    S.MixingParameterSimpleMag = S.MixingParameterMag;
end
    
% Preconditioner for SCF convergence
if S.MixingVariable < 0
    S.MixingVariable = 0; % set default mixing var to density
end

if S.MixingPrecond < 0
    S.MixingPrecond = 1; % set default precond to 'Kerker' preconditioner
end

if S.MixingPrecondMag < 0
    S.MixingPrecondMag = 0; % set default precond to none
end

% set up coefficients
if S.MixingPrecond == 1 % kerker
    S.precondcoeff_a = 1.0;
    S.precondcoeff_lambda_TF = S.precond_kerker_kTF * S.precond_kerker_kTF;
    S.precondcoeff_k = 0;
elseif S.MixingPrecond == 2 % resta
    % put these in input options
    %max_q = sqrt(3) * 2*pi/0.08;  % Maximum q in the domain of fit
    max_q = 100;
    % mpower = 2; % The power used in the fit 
    % a0 = 0.25;   % parameters of the preconditioner
    % ktf = 1;     % parameters of the preconditioner
    % q0 = 1.36;  % parameters of the preconditioner
    % e0 = 5.7;   % parameters of the preconditioner
    % Rs = 2.76;  % parameters of the preconditioner
    mpower = S.precond_fitpow;
    ktf    = S.precond_kerker_kTF;
    a0     = S.precond_kerker_thresh;
    q0     = S.precond_resta_q0;
    Rs     = S.precond_resta_Rs;
    e0 = sinh(q0*Rs)/(q0*Rs);   % parameters of the preconditioner
    %[a_temp, lambda_TF_temp, const_temp] = fit_mixing_preconditioner(...
    %    100,1.36,5.7,2.76,0.25,1,2,2);
    [a_temp, lambda_TF_temp, const_temp] = fit_mixing_preconditioner(...
        max_q, q0, e0, Rs, a0, ktf, mpower, 2);
    % store coeffs in S
    S.precondcoeff_a         = a_temp;
    S.precondcoeff_lambda_TF = lambda_TF_temp;
    S.precondcoeff_k         = const_temp;
elseif S.MixingPrecond == 3 % truncated kerker
    % put these in input options
    %max_q = sqrt(3) * 2*pi/0.08;  % Maximum q in the domain of fit
    max_q = 100;
    %mpower = 2; % The power used in the fit 
    %a0 = 0.25;   % parameters of the preconditioner
    %ktf = 1;     % parameters of the preconditioner
    %q0 = 1.36;   % parameters of the preconditioner
    %Rs = 2.76;   % parameters of the preconditioner
    mpower = S.precond_fitpow;
    ktf    = S.precond_kerker_kTF;
    a0     = S.precond_kerker_thresh;
    q0     = S.precond_resta_q0;
    Rs     = S.precond_resta_Rs;
    e0 = sinh(q0*Rs)/(q0*Rs);   % parameters of the preconditioner
    [a_temp, lambda_TF_temp, const_temp] = fit_mixing_preconditioner(...
        max_q, q0, e0, Rs, a0, ktf, mpower, 1);
    S.precondcoeff_a         = a_temp;
    S.precondcoeff_lambda_TF = lambda_TF_temp;
    S.precondcoeff_k         = const_temp;
end

if (S.RelaxFlag || S.MDFlag)
    % Name of the restart file
    S.restartfname = strcat(filename,'.restart');
    % charge difference histories for charge extrapolation
    S.delta_rho_tm2 = zeros(S.N,1);
    S.delta_rho_tm1 = zeros(S.N,1);
    S.delta_rho_t   = zeros(S.N,1);
    S.delta_rho_in_tp1 = zeros(S.N,1);
    % S.dV_tm2 = S.dV;
    % S.dV_tm1 = S.dV;
    % S.dV_t   = S.dV;
    % S.dV_tp1 = S.dV;
    S.atom_pos_tm2  = S.Atoms_init;
    S.atom_pos_tm1  = S.Atoms_init;
    S.atom_pos_t    = S.Atoms_init;
    S.atom_pos_tp1  = S.Atoms_init;
    S.Atoms_old     = S.Atoms_init;
end

if S.rhoTrigger < 0
    if S.spin_typ == 2
        S.rhoTrigger = 6;
    else
        S.rhoTrigger = 4;
    end
end

S.Relax_iter = 1;
S.ForceCount = 1;
S.amu2au = 1822.888485332371; % 1 au = 9.10938356e-31 Kg; 1 amu =  1.660539040e-27 Kg;
S.fs2atu = 41.34137333649300; %1atu = 2.418884326509e-17 s;

fprintf(' Creating differentiation matrices ...\n');
t1 = tic;

% Calculate discrete laplacian (1D) and discrete gradient indices' values
S = lapIndicesValues_1d(S);
S = gradIndicesValues(S);

% Calculate discrete laplacian
[DL11,DL22,DL33,DG1,DG2,DG3] = blochLaplacian_1d(S,[0 0 0]);
if S.cell_typ < 3
    S.Lap_std = S.lapc_T(1,1) * kron(speye(S.Nz),kron(speye(S.Ny),DL11))  +  S.lapc_T(2,2) * kron(speye(S.Nz),kron(DL22,speye(S.Nx))) + ...
                S.lapc_T(3,3) * kron(DL33,kron(speye(S.Ny),speye(S.Nx))) ;
    if (S.cell_typ == 2)
        MDL = S.lapc_T(1,2) * kron(speye(S.Nz),kron(DG2,DG1)) + S.lapc_T(2,3) * kron(DG3,kron(DG2,speye(S.Nx))) + ...
              S.lapc_T(1,3) * kron(DG3,kron(speye(S.Ny),DG1)) ;
        S.Lap_std = S.Lap_std + MDL;
    end
elseif (S.cell_typ == 3 || S.cell_typ == 4 || S.cell_typ == 5)
    S.Lap_std = kron(speye(S.Nz),kron(speye(S.Ny),(DL11+DG1))) + kron(DL33,kron(speye(S.Ny),speye(S.Nx))) + ...
                kron(speye(S.Nz),kron(DL22,S.R2inv));
    if (S.cell_typ == 4 || S.cell_typ == 5)
        MDL = kron(speye(S.Nz),kron(DL22,speye(S.Nx))) +   kron(DG3,kron(DG2,speye(S.Nx)));
        S.Lap_std = S.Lap_std + MDL;
    end
end

% Calculate discrete gradient
S.grad_1 = blochGradient(S,[0 0 0],1);
S.grad_2 = blochGradient(S,[0 0 0],2);    
S.grad_3 = blochGradient(S,[0 0 0],3);

% Calculate preconditioners for negative discrete laplacian
[S.LapPreconL, S.LapPreconU] = ilu(S.Lap_std,struct('droptol',1e-5));

% initialize vdWDF
if (S.vdWDFFlag == 1) || (S.vdWDFFlag == 2) % 1: temporary flag of vdW-DF1 2: vdW-DF2
    if S.BC ~= 2 % vdWDF can only be used in 3D periodic boundary condition because it used FFT
        error('vdW-DF can only be used in 3D periodic system!');
    end
        
    if S.vdWDFKernelGenFlag == 1
        S = vdWDF_Initial_GenKernel(S);
    else % input the saved Kernel function for saving time
        S = vdWDFinitialize_InputKernel(S);
    end
end

fprintf(' Done. (%.3f sec)\n', toc(t1));

% Estimate memory usage
S.memory_usage = my_estimate_memory(S);

if S.usefock == 1
    S = my_exx_initialization(S);
end

end