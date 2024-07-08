function S = my_inpt_defaults()
% Set up the default parameters

% cell_typ: 1 - orthogonal unit cell or 2 - non-orthogonal unit cell
cell_typ = 1; % default is orthogonal

% lat_vec: Lattice unit vectors
lat_vec = eye(3);
% corresponding metric_T, grad_T, lapc_T and Jacb
metric_T = eye(3);
grad_T   = eye(3);
lapc_T   = eye(3);
Jacb     = 1;

% L1, L2, L3: domain size in each direction
% no defaults
L1 = 0.0;
L2 = 0.0;
L3 = 0.0;

% Nx, Ny, Nz: number of finite-difference grids in each direction
% no defaults
Nx = 0;
Ny = 0;
Nz = 0;
% N: N = Nx * Ny * Nz
N  = 0;

% Ecut
ecut = -1.0;

% mesh spacing
mesh_spacing = -1.0;

% kptgrid: k-point grid
kptgrid = [0.0 0.0 0.0]; % default is Gamma-point

% shift in k-point grid    
kptshift = [0.0 0.0 0.0];

% nkpt: number of k-points
nkpt = [1 1 1];

% tnkpt: number of Time-Reversal Symmetry reduced k-points
tnkpt = 1;

% wkpt: weights for k-points
wkpt = 1;

% BC: boundary conditions. 
% 1 - Zero Dirichlet, 2 - Periodic, 
% 3 - surface in x-y, 4 - wire along z-direction
BC = -1; % will be set up later after reading user inputs
BCx = -1; BCy = -1; BCz = -1;

% isBS: flag for band structure calculation
isBS = 0; % Default is off

% lattice: lattice type, for band structure calculation
lattice = 'undefined';

% SCF_tol: SCF tolerance
SCF_tol = -1.0; % default will be set after inpt is read

% Nev: Number of states/bands
Nev = -1; % default will be set after inpt is read

% poisson_tol: Poisson tolerance
poisson_tol = -1; % default will be set after inpt is read

% pseudocharge_tol: Pseudocharge (rb) tolerance
pseudocharge_tol = -1;

% Cst: Factor for conversion from Ha to eV
Cst = 27.21138602;

% Temp: Electronic temperature, beta := 1/ (kB * Temp)
%Temp = 315.7751307269723;
kB = (8.6173303e-5)/Cst;
%bet = 1.0 / (kB * Temp); % beta = 1 / smearing
elec_T_type = 1; % gaussian smearing
bet = -1;
Temp = -1;

% npl: Degree of Chebyshev polynomial
npl = -1; % default will be set after inpt is read

% conditioner: mixing type and preconditioner
% 1 - Potential mixing , 2 - Density + Kerker, 3 - Density mixing
% conditioner = 1;

% FDn: half of finite difference order
FDn = 6;

% rc_ref: rc reference
rc_ref = 0.5;

% max_relax_it: Maximun number of relaxations
max_relax_it = 0;

% dbg_switch: debug switch
dbg_switch = 0;

% xc: Exchange-correlation functional 
% 0 - LDA_PW, 1 - LDA_PZ, 2 - PBE(GGA)
xc = 0;

% Nelectron: number of electrons
% no default
Nelectron = 0;

% n_typ: number of atom types
n_typ = 0;

S = struct(...
    'cell_typ',cell_typ,'lat_vec',lat_vec,'metric_T',metric_T,...
    'grad_T',grad_T, 'lapc_T',lapc_T,'Jacb',Jacb,'L1',L1,'L2',L2,'L3',L3,...
    'Nx',Nx,'Ny',Ny,'Nz',Nz,'N',N,'ecut',ecut,'mesh_spacing',mesh_spacing,'kptgrid',kptgrid,...
    'kptshift',kptshift,'nkpt',nkpt,'tnkpt',tnkpt,'wkpt',wkpt,'BC',BC,'BCx',BCx,'BCy',BCy,'BCz',BCz,...
    'isBS',isBS,'lattice',lattice,'SCF_tol',SCF_tol,'Nev',Nev,'poisson_tol',poisson_tol,...
    'pseudocharge_tol',pseudocharge_tol, 'Cst',Cst,'kB',kB,'elec_T_type',elec_T_type,...
    'Temp',Temp,'bet',bet,'npl',npl,'FDn',FDn,...
    'rc_ref',rc_ref,'max_relax_it',max_relax_it,...
    'dbg_switch',dbg_switch,'xc',xc,...
    'Nelectron',Nelectron,'n_typ',n_typ);

S.TimeRevSym = 1; 

S.spin_typ = 0;
S.nspin = 1; % spin 

% EXTRA variables from SPARC
S.CheFSI_Optmz = 0;
S.chefsibound_flag = 0;
S.FixRandSeed = 0;
S.rhoTrigger = -1;
S.nchefsi = 1;
S.NetCharge = 0;
S.MAXIT_SCF = 100;
S.MINIT_SCF = 2;
S.MAXIT_POISSON = 1000;
S.accuracy_level = -1;
S.target_force_accuracy = -1.0;
S.target_energy_accuracy = -1.0;
S.TOL_RELAX = 5e-4;
S.TOL_LANCZOS = 1e-2;
S.StandardEigenFlag = 0;

% preconditioning
S.precond_tol = -1;

S.precond_kerker_kTF = 1;
S.precond_kerker_thresh = 0.1;
S.precond_kerker_kTF_mag = 1;
S.precond_kerker_thresh_mag = 0.1;
% S.precond_fitpow = 2;
% S.precond_resta_q0 = 1.36;
% S.precond_resta_Rs = 2.76;

S.MixingVariable = -1;
S.MixingPrecond = -1;
S.MixingPrecondMag = -1;
S.Pf_guess = [];

S.MixingHistory = 7;
S.MixingParameter = 0.3;
S.MixingParameterSimple = -1.0; % for simple mixing, set up later
S.MixingParameterMag = -1.0; % default mixing parameter for magnetization density/potential
S.MixingParameterSimpleMag = -1.0; % for simple mixing, set up later
S.PulayFrequency = 1;
S.PulayRestartFlag = 0;

S.TWtime = 1000000000;
S.RelaxFlag = 0;
S.RelaxMeth = 'LBFGS';
S.max_relax_it = 100;
S.max_dilatation = 1.2;    
S.TOL_RELAX_CELL = 0.01; % in GPa (max pressure)
S.MDFlag = 0;
S.RestartFlag = 0;
S.MDMeth = 'NVE';
S.MD_dt = 1.0;
S.MD_Nstep = 0;
S.ion_T = -1.0;
S.thermos_TF = -1.0;
S.ion_elec_eqT = 1;
S.ion_vel_dstr = 2;
S.NLCG_sigma = 0.5;
S.qmass = 1;
S.L_history = 20;
S.L_finit_stp = 5e-3;
S.L_maxmov = 0.2;
S.L_autoscale = 1;
S.L_lineopt = 1;
S.L_icurv = 1.0;
S.FIRE_dt = 1.0;
S.FIRE_mass = 1.0;
S.FIRE_maxmov = 0.2;
S.Calc_stress = 0;
S.Calc_pres = 0;
S.PrintForceFlag = 1;         % flag for printing forces
S.PrintAtomPosFlag = 1;       % flag for printing atomic positions
S.PrintAtomVelFlag = 1;       % flag for printing atomic velocities
S.PrintElecDensFlag = 0;      % flag for printing final electron density
S.PrintEigenFlag = 0;         % Flag for printing final eigenvalues and occupations
S.PrintMDout = 1;             % Flag for printing MD output in a .aimd file
S.PrintRelaxout = 1;          % Flag for printing relax output in a .relax file
S.Printrestart = 1;           % Flag for printing output needed for restarting a simulation
S.Printrestart_fq = 1;        % Steps after which the output is written in the restart file

S.vdWDFFlag = 0;              % Flag for calculating vdW-DF
S.vdWDFKernelGenFlag = 0;     % Flag for calculating kernel functions of vdW-DF

% Cell option
S.Flag_latvec_scale = 0;
S.latvec_scale_x = 0.0;
S.latvec_scale_y = 0.0;
S.latvec_scale_z = 0.0;

% NLCC
S.NLCC_flag = 0;

% Origin of the unit cell wrt some global origin
S.xin = 0;  
S.yin = 0;
S.zin = 0;

% Cychel
S.alph = 0.0;

% SOC
S.SOC_flag = 0;
S.nspinor = 1;
S.nspden = 1;

% DFT-D3 parameters
S.d3Flag = 0;
S.d3Rthr = 1600.0;
S.d3Cn_thr = 625.0;

% hybrid functionals
S.usefock = 0;
S.MAXIT_FOCK = -1;
S.MINIT_FOCK = -1;
S.FOCK_TOL = -1;
S.hyb_mixing = 0.0;
S.hyb_range_fock = -1;
S.hyb_range_pbe = -1;
S.ExxMethod = '';
S.SCF_tol_init = -1;
S.ACEFlag = 1;
S.EXXACEVal_state = 3;
S.exx_downsampling = [1 1 1];
S.ExxDivMethod = '';
end