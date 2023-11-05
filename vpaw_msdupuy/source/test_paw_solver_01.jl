Pkg.activate("../")
push!(LOAD_PATH, pwd())

import Eigensolvers

import PWCoulombModule
using PWCoulombModule: PWSolverCoulomb
const pw_coulomb = PWCoulombModule
const params = PWSolverCoulomb

import VPAWVanderbiltModule
const paw = VPAWVanderbiltModule

import VPAWSolverModule
const vpawsolver = VPAWSolverModule

import PAWSolverModule
const paw_solver = PAWSolverModule

function test_fft_H(N,L,rc,Npaw,Z)
   coef_TM = paw.coef_TM(rc, 1, 0, Z, 1e-8)[1]
   V(x,y,z) = paw.V_scr(norm([x,y,z]-[L/2,L/2,L/2]), 1, 0, rc, coef_TM,Z)
   p = params(N,L,V)
   X = zeros(3,1)
   X[:,1] = [L/2, L/2, L/2]
   fpaw = pawfunc(rc, X, p, Npaw, Z)
   psi, E, res = energy_paw(fpaw, p, vpawsolver.tdphi_test(N,L,[L/2,L/2,L/2], p, rc,Z))
   return psi, E
end

function test_num_H(N,L,rc,Npaw,Z;args...)
   coef_TM = paw.coef_TM(rc, 1, 0, Z, 1e-8)[1]
   V(r) = paw.V_scr(r, 1, 0, rc, coef_TM, Z) 
   p = params(N,L,X,Z,V)
   X = zeros(3,1)
   X[:,1] = [L/2, L/2, L/2]
   fpaw = pawfunc(rc, X, p, Npaw, Z; proj=vpawsolver.proj_num, args...)
   psi, E, res = energy_paw(fpaw, p, vpawsolver.tdphi_test(N,L,[L/2,L/2,L/2], p, rc,Z))
   return psi, E
end

function test_fft_H2(N,L,rc,Npaw,R,Z)
   X1 = [(L-R)/2,L/2,L/2]
   X2 = [(L+R)/2,L/2,L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   coef_TM = paw.coef_TM(rc, 1, 0, Z, 1e-8)[1]
   V(x,y,z) = paw.V_scr(norm([x-(L-R)/2,y-L/2,z-L/2]), 1, 0, rc, coef_TM, Z) + paw.V_scr(norm([x-(L+R)/2,y-L/2,z-L/2]), 1, 0, rc, coef_TM, Z)
   p = params(N,L,V)
   fpaw = pawfunc(rc, X, p, Npaw, Z)
   psi, E, res = energy_paw(fpaw, p, vpawsolver.guess_H2(rc,X1,X2,p))
   return psi, E
end

function test_num_H2(N,L,rc,Npaw,R,Z;args...)
   X1 = [(L-R)/2,L/2,L/2]
   X2 = [(L+R)/2,L/2,L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   coef_TM = paw.coef_TM(rc, 1, 0, Z, 1e-8)[1]
   V(r) = paw.V_scr(r, 1, 0, rc, coef_TM, Z);
   p = params(N,L,X,Z,V);
   fpaw = paw_solver.pawfunc(rc, X, p, Npaw, Z; proj=vpawsolver.proj_num, args...);
   psi, E, res = paw_solver.energy_paw(fpaw, p, vpawsolver.guess_H2(rc,X1,X2,p,Z))
   return psi, E
end

N = 10
L = 10.0
rc = 2.0 # cutoff radius?
Npaw = [1, 1]
R = 1.2 # bond length
Z = 1.0
psi, E = test_num_H2(N,L,rc,Npaw,R,Z)
