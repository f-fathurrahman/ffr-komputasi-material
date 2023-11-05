module PAWSolverModule

using Test
using QuadGK
using Plots
using Polynomials.PolyCompat
using GSL
using LinearAlgebra
using FFTW

import Eigensolvers

import PWCoulombModule
using PWCoulombModule: PWSolverCoulomb
const pw_coulomb = PWCoulombModule
const params = PWSolverCoulomb

import VPAWVanderbiltModule
const paw = VPAWVanderbiltModule

import VPAWSolverModule
const vpawsolver = VPAWSolverModule


mutable struct pawfunc
   rc :: Float64 #cut-off radius
   p :: params
   X :: Array{Float64, 2} #positions of the atoms
   Nat :: Integer #number of atoms
   Npaw :: Array{Int64,1} #number of PAW functions
   Npawtot :: Integer
   P #projectors
   DH #Array of matrices <phi - \tilde\phi | H  |phi - \tilde\phi >
   DS #Array of matrices <phi - \tilde\phi | phi - \tilde\phi >
   function pawfunc(rc, X::Array{Float64,2}, p::params, Npaw:: Array{Int64,1}, Z; proj=vpawsolver.proj_fft, args...)
      fpaw = new(rc, p, X, size(X)[2], Npaw, paw.Npawtot(Npaw), nothing)
      Npawtot = paw.Npawtot(Npaw)
      coefpaw = paw.pawcoef(Z, rc, Npaw; args...)
      fpaw.DH = zeros(Float64, (Npawtot,Npawtot,fpaw.Nat))
      fpaw.DS = zeros(Float64, (Npawtot,Npawtot,fpaw.Nat))
      fpaw.P = zeros(ComplexF64, (p.size_psi..., Npawtot, fpaw.Nat))
      for iat in 1:fpaw.Nat
         fpaw.DS[:,:,iat] = S_ij(rc, Z, coefpaw.tdR, Npaw)
         fpaw.DH[:,:,iat] = D_ij(rc, X, coefpaw.tdR, iat, coefpaw.coef_TM, Npaw, Z)[1]
         fpaw.P[:,:,:,:,iat] = proj(rc, X[:,iat], p::params, coefpaw::paw.pawcoef)
      end
      return fpaw
   end
end

#-\Delta(R_nl R_nl - tdR_nl td_Rnl
function Ddiff(r, rc, n1, n2, l, coef_PAW,  Z)
   if r > rc
      return 0.
   else
      L(r,n) = GSL.sf_laguerre_n(n-l-1, 2l+1, 2Z*r/n)
      function lag_der(r,n) #derivative of the associated Laguerre polynomials
         if 1 > n-l-1
            return 0.
         else
            return (-1)*GSL.sf_laguerre_n(n-l-2, 2+2l, 2Z*r/n)
         end
      end
      der_phi(r,n) = paw.C_nl(n,l,Z)*exp(-Z*r/n)*(-Z*r^(l+1)*L(r,n)/n + (l+1)*r^l*L(r,n) + 2Z/n*lag_der(r,n)*r^(l+1))
      der_tdphi(r,n) = paw.C_nl(n,l,Z)*polyval(polyder(Poly(coef_PAW[:,n-l,l+1])(Poly([0,0,1]))*poly(zeros(l+1))), r)
      return der_phi(r,n1)*der_phi(r,n2) - der_tdphi(r,n1)*der_tdphi(r,n2)
   end
end

function D_ij(rc, X, coef_PAW, N, coef_TM, Npaw, Z) #N : atom site (working)
   Rnl_diff(r, rc, l1, l2, coef_PAW, j, k) = paw.R_nl(r, l1+j, l1, Z)*paw.R_nl(r, l2+k, l2, Z) - paw.tilde_R_nl(r, rc, l1+j, l1, Z, coef_PAW[:,j,l1+1])*paw.tilde_R_nl(r, rc, l2+k, l2, Z, coef_PAW[:,k,l2+1])
   Npawtot = paw.Npawtot(Npaw)
   D = zeros(ComplexF64,(Npawtot,Npawtot))
   err = zeros(ComplexF64,(Npawtot,Npawtot))
   ind = 0
   for lpaw in eachindex(Npaw)
      l = lpaw -1
      for m in 1:(2lpaw-1)
         for j in 1:Npaw[lpaw]
            for k in 1:Npaw[lpaw]
               D[j+ind,k+ind] = 0.5*QuadGK.quadgk(r -> Ddiff(r, rc, j+l, k+l, l, coef_PAW, Z), 0, rc)[1]
               D[j+ind,k+ind] += 0.5*QuadGK.quadgk(r -> l*(l+1)/r^2*Rnl_diff(r, rc, l, l, coef_PAW, j, k), 0, rc)[1]
               D[j+ind,k+ind] -= QuadGK.quadgk(r -> Z/r*paw.R_nl(r, l+j, l, Z)*paw.R_nl(r, l+k, l, Z), 0, rc)[1] #Z/r \phi_nl \phi_n'l
               D[j+ind,k+ind] -= QuadGK.quadgk(r -> paw.V_scr(r,1,0,rc,coef_TM,Z)*paw.tilde_R_nl(r, rc, l+j, l, Z, coef_PAW[:,j,l+1])*paw.tilde_R_nl(r, rc, l+k, l, Z, coef_PAW[:,k,l+1]), 0, rc)[1]
               #V_TM \tilde \phi_nl \tilde\phi_n'l
            end
         end
         ind += Npaw[lpaw]
      end
   end
   if size(X)[2]==1
      return D, err
   else
      Y = X[:,1:end .!=N]
      pot(x,y,z) = Z*sum(-1/norm([x,y,z] + X[:,N] - Y[:,i3]) for i3 in 1:(length(X[1,:])-1))
      for i in 1:Npawtot
         for j in 1:Npawtot
            l1, n1, m1 = vpawsolver.int_to_nl(i, Npaw)
            l2, n2, m2 = vpawsolver.int_to_nl(j, Npaw)
            if l1 == l2 == 1 #use radial symmetry to recast to 1D integral
               for i3 in 1:(size(X)[2]-1)
                  R =  norm(X[:,N]-Y[:,i3])
                  pot2, err2 = QuadGK.quadgk(r -> Z/R*Rnl_diff(r, rc, l1-1, l2-1, coef_PAW, n1, n2), 0, rc)
                  D[i,j] += - pot2
                  err[i,j] += err2
               end
            else  #other l1 and l2 values
               realintg(X3) = real( pot(X3...)*Rnl_diff(norm(X3), rc, l1-1, l2-1, coef_PAW, n1, n2)*conj(paw.Y_lm(X3...,l1-1,m1))*paw.Y_lm(X3...,l2-1,m2))/norm(X3)^2
               imagintg(X3) = imag( pot(X3...)*Rnl_diff(norm(X3), rc, l1-1, l2-1, coef_PAW, n1, n2)*conj(paw.Y_lm(X3...,l1-1,m1))*paw.Y_lm(X3...,l2-1,m2))/norm(X3)^2
               rea, rea_err = vpawsolver.int3D(realintg,rc)
               D[i,j] += rea
               err[i,j] += rea_err
               ima, im_err = vpawsolver.int3D(imagintg,rc)
               D[i,j] += im*ima
               err[i,j] += im*im_err
            end
         end
      end
      return D, err
   end
end

function S_ij(rc, Z, coef_PAW, Npaw)
   S = zeros(paw.Npawtot(Npaw), paw.Npawtot(Npaw))
   ind = 0
   Rnl_tdRnl(r,lpaw,j,k) = paw.R_nl(r, j+lpaw-1, lpaw-1, Z)*paw.R_nl(r, k+lpaw-1, lpaw-1, Z) - paw.tilde_R_nl(r, rc, j+lpaw-1, lpaw-1, Z, coef_PAW[:,j,lpaw])*paw.tilde_R_nl(r, rc, k+lpaw-1, lpaw-1, Z, coef_PAW[:,k,lpaw])
   for lpaw in eachindex(Npaw)
      for m in 1:(2lpaw-1)
         for j in 1:Npaw[lpaw]
            for k in 1:Npaw[lpaw]
               S[j+ind,k+ind] = QuadGK.quadgk(r -> Rnl_tdRnl(r,lpaw,j,k), 0., rc, atol=1e-10, rtol=1e-10)[1]
            end
         end
         ind += Npaw[lpaw]
      end
   end
   return S
end

"""
tdH_paw computes tilde H tilde psi
"""

function tdH_paw(fpaw::pawfunc, p::params, psi)
   L = p.L1*p.L2*p.L3
   Ppsi = zeros(ComplexF64, (fpaw.Npawtot, fpaw.Nat))
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         Ppsi[ipaw, iat] = 1/p.Ntot^2*L*dot(fpaw.P[:,:,:,ipaw,iat], psi)
      end
   end
   tdpsi = PWSolverCoulomb.ham(p, psi)
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         tdpsi +=  (fpaw.DH[:,:,iat]*Ppsi[:,iat])[ipaw]*fpaw.P[:,:,:,ipaw,iat]
      end
   end
   return tdpsi
end

"""
tdH_paw computes tilde S tilde psi
"""

function tdS_paw(fpaw::pawfunc, p::params, psi)
   L = p.L1*p.L2*p.L3
   Ppsi = zeros(ComplexF64, (fpaw.Npawtot, fpaw.Nat))
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         Ppsi[ipaw,iat] = 1/p.Ntot^2*L*dot(fpaw.P[:,:,:,ipaw,iat], psi)
      end
   end
   tdpsi = psi
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         tdpsi += (fpaw.DS[:,:,iat]*Ppsi[:,iat])[ipaw]*fpaw.P[:,:,:,ipaw,iat]
      end
   end
   return tdpsi
end

function energy_paw(fpaw::pawfunc, p::params, seed)
   H_1var(psi) = tdH_paw(fpaw, p, reshape(psi, p.size_psi))[:]
   S_1var(psi) = tdS_paw(fpaw, p, reshape(psi, p.size_psi))[:]
   function P(psi)
       meankin = sum(p.kin[i]*abs2(psi[i]) for i = 1:p.Ntot) / (norm(psi)^2)
       return psi ./ (0.1*meankin .+ p.kin[:]) # this should be tuned but works more or less
    end
   return Eigensolvers.eig_lanczos(H_1var, seed[:], B=S_1var, m=3, Imax = 250, do_so=true, norm_A = 6pi^2*p.N1)
   #return eigensolvers.eig_pcg(H_1var, seed[:];P =P, B=S_1var, tol=1e-3, maxiter = 1000)
end


end # module
