% This is truncated from electronicGroundStateAtomicForce

% Must call setup_* first, i.e. S struct must be available

% Set an argument
atom_pos = S.Atoms;


%
% Here comes the actual implementation
%

assert(S.RelaxFlag == 0, 'RelaxFlag must be 0 for this script')
assert(S.MDFlag == 0, 'MDFlag must be 0 for this script')

% Reference_cutoff checks
nndis = my_calculate_min_distance(S);
if S.rc_ref > 0.5*nndis
    fprintf("WARNING: neighbor distance is too small\n")
end
if S.rc_ref < S.dx || S.rc_ref < S.dy || S.rc_ref < S.dz
    fprintf("WARNING: Too large mesh spacing\n")
end
    
% Check position of atom near the boundary and apply wraparound in case of PBC
S = my_check_atomlocation(S);

% Pseudocharge (and reference), sum atomic charge density, self energy 
% (and reference), electrostatic correction 
% S.b,S.b_ref,S.Eself,S.Eself_ref,S.rho_at,S.E_corr,S.V_c, S.NegCharge,
% S.PosCharge, S.NetCharge
S = calculate_b_guessRho_Eself(S);

% set up guess electron density (guess rho)
S = my_initElectrondensity(S);

% Calculate nonlocal projectors    
S.Atom = calculate_nloc_projector(S);

% Self-consistent Field (SCF) method
S = my_scf(S);

