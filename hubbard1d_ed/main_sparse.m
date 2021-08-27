clear all;
close all;

%% For running on SciClone
format compact;

Lx = 2; % The number of lattice sites in the x direction
N_up = 1; % The number of spin-up electrons
N_dn = 1; % The number of spin-down electrons
U = 1.1; % The on-site repulsion strength in the Hubbard Hamiltonian
t = 1.2; % The hopping amplitude between nearest-neighbor sites in the x direction

tau = 0.0;

maxNumCompThreads(2)

[ exactUp, exactDn ] = unequalTimeGF( t, U, tau, Lx, N_up, N_dn );

exactUp

exactDn
