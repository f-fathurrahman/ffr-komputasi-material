function memory_usage = estimate_memory(S)
% estimate the memory required for the simulation
size_double = 8;

fd_order = S.FDn * 2;

N_ex = (S.Nx+fd_order) * (S.Ny+fd_order) * (S.Nz+fd_order);

% orbitals (dominant)
ncpy_orbitals = 3; % 4 copies required during chebyshev filtering
if S.nspin ~= 1, ncpy_orbitals = ncpy_orbitals * 2; end
% for kpoints, the factor 2 is for complex entries
if S.tnkpt ~= 1, ncpy_orbitals = ncpy_orbitals * 2 * S.tnkpt; end
if S.nspinor ~= 1, ncpy_orbitals = ncpy_orbitals * 2; end
memory_orbitals = S.N * S.Nev * size_double * ncpy_orbitals;

% sparse matrices
% Lap_std, grad_i, {GIi, GJi, GVi}
% (nnz * 2 + ncol + 1) * size_double
nnz_sparse = (2 * (3*fd_order+1) + 3 * fd_order) * S.N;
memory_sparse = (2*nnz_sparse + S.N + 1) * size_double ...
    + 9*fd_order*S.N * size_double + (S.BC ~= 1) * 6*fd_order*S.N;

% vectors, rho, phi, Veff, ...; isIn, RR_AUG, RR_AUG_3D
memory_vectors = 14 * S.N * size_double + ...
    + (S.BC == 1) * 3 * N_ex * size_double ...
    + (S.RelaxFlag || S.MDFlag) * 4 * S.N * size_double;
    
% spherical harmonics
if S.BC == 1
    memory_SH = (6+1)^2 * (S.N + N_ex) * size_double;
else 
    memory_SH = 0;
end

% history matrix
memory_hist = 2 * S.MixingHistory * S.N * size_double;

% total
memory_usage = memory_orbitals + memory_sparse + memory_vectors ...
    + memory_SH + memory_hist;

fprintf('\n');
fprintf(' Estimated memory usage:\n');
fprintf(' Total: %s\n', my_print_mem(memory_usage));
fprintf(' orbitals            : %s\n', my_print_mem(memory_orbitals));
fprintf(' sparse matrices     : %s\n', my_print_mem(memory_sparse));
fprintf(' global-size vectors : %s\n', my_print_mem(memory_vectors));
if S.BC == 1
    fprintf(' spherical harmonics : %s\n', my_print_mem(memory_SH));
end
fprintf(' mixing histories    : %s\n', my_print_mem(memory_hist));
fprintf('\n');

end
