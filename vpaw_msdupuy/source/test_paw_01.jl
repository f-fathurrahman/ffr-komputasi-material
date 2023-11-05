Pkg.activate("../")
push!(LOAD_PATH, pwd())

using Test
using Plots
using Polynomials
using LinearAlgebra

import Eigensolvers

import PWCoulombModule
using PWCoulombModule: PWSolverCoulomb
const pw_coulomb = PWCoulombModule
const params = PWSolverCoulomb

import VPAWVanderbiltModule
const paw = VPAWVanderbiltModule

import VPAWSolverModule
const vpawsolver = VPAWSolverModule



"""
Orthogonality of  tilde phi and  tilde p

tdphi_fft to test orthogonality of  tilde phi_j and  tilde p_k
"""
function ortho_test(rc,N,L,Npaw,Z,mult;proj=vpawsolver.proj_num)
   V(x,y,z) = 1.
   X = zeros(3,1)
   p = params(N,L,X,Z)
   coefpaw = paw.pawcoef(Z, rc, Npaw,GS = paw.GS_custom, proj_gen = paw.coef_proj_custom)
   X = zeros(3,1)
   X[:,1] = [L/2, L/2, L/2]
   P = proj(rc, X[:,1], p, coefpaw)
   Npawtot = paw.Npawtot(Npaw)
   tdphi(r,l,n) = paw.tilde_R_nl(r,rc,n+l-1,l-1,Z,coefpaw.tdR[:,n,l])
   temp_tdphi = zeros(ComplexF64, (2p.N1*mult+1, 2p.N2*mult+1, 2p.N3*mult+1, Npawtot))
   tdPhi = zeros(ComplexF64, (p.size_psi..., Npawtot))
   for ipaw in 1:(Npawtot)
      l,n,m = vpawsolver.int_to_nl(ipaw,Npaw)
      for i1 in 1:(2*p.N1*mult+1)
         for i2 in 1:(2*p.N2*mult+1)
            for i3 in 1:(2*p.N3*mult+1)
               #r = norm(vpawsolver.coords(p,i1,i2,i3,mult) - X[:,1])
               r = sqrt( ((i1-1)/(2*p.N1*mult+1)*p.L1 - L/2)^2 + ((i2-1)/(2*p.N2*mult+1)*p.L2 -L/2)^2 + ((i3-1)/(2*p.N3*mult+1)*p.L3 - L/2)^2 )
               temp_tdphi[i1,i2,i3,ipaw] = paw.Y_lm(pw_coulomb.coords(p,i1,i2,i3,mult) - X[:,1]..., l-1, m)*tdphi(r,l,n)/r
            end
         end
      end
      @views temp_tdphi[:,:,:,ipaw] = fft(temp_tdphi[:,:,:,ipaw]) #cannot do on site fft (not doing the right thing)
      @views tdPhi[:,:,:,ipaw] = pw_coulomb.fft_reshape(temp_tdphi[:,:,:,ipaw], p, mult, 1)
   end
   out = zeros(ComplexF64,(Npawtot,Npawtot))
   for i1 in 1:Npawtot
      for i2 in 1:Npawtot
         out[i1,i2] = 1/p.Ntot^2*L^3*vecdot(P[:,:,:,i1],tdPhi[:,:,:,i2])
      end
   end
   return out
end

#return \psi = (I+T) \tilde\psi : TO BE WRITTEN
function tdpsi_to_psi(tdpsi,rc,N,L,Npaw,Z,mult;proj=vpawsolver.proj_num) #tdpsi in reciprocal space
   coefpaw = paw.pawcoef(Z, rc, Npaw)
   X = zeros(3,1)
   X[:,1] = [0.5L,0.5L,0.5L]
   p = params(N,L,X,Z)
   P = proj(rc,X[:,1],p,coefpaw)
end

