
module VPAWSolverModule

using FFTW
using Test
using QuadGK
using Plots
using Polynomials
using LinearAlgebra
using Combinatorics
using Cubature
using GSL

import Eigensolvers

import PWCoulombModule
using PWCoulombModule: PWSolverCoulomb
const pw_coulomb = PWCoulombModule
const params = PWSolverCoulomb

import VPAWVanderbiltModule
const paw = VPAWVanderbiltModule


"""
td : means tilde
"""

function fft_mode(i,N)
   if (i <= N+1)
      return i-1
   else
      return i-2*N-2
   end
end

#converts integer to (n,l,m) quantum numbers for PAW functions
#compared to usual quantum number, int_to_nl returns : [l+1,n-l-1,m]
function int_to_nl(ipaw,Npaw)
   l = 0
   while (ipaw > (2l+1)*Npaw[l+1])
      ipaw -= (2l+1)*Npaw[l+1]
      l += 1
   end
   m = -l
   while (ipaw/Npaw[l+1] > 1)
      ipaw -= Npaw[l+1]
      m += 1
   end
   n = ipaw
   return [l+1,n,m]
end

function int3D(f,rc) #integration on the box [-rc,rc]^3 for functions with singularity at (0,0,0)
   out = 0.
   err = 0.
   for i1 in 0:1
      for i2 in 0:1
         for i3 in 0:1
            x,y = Cubature.hcubature(f, [-rc*(1-i1),-rc*(1-i2),-rc*(1-i3)], [rc*i1,rc*i2,rc*i3],abstol=1e-7,reltol=1e-7,maxevals=100000)
            out += x
            err += y
         end
      end
   end
   return out, err
end

"""
create the matrix P of size M(==(2N1+1)(2N2+1)(2N3+1)) x N_proj
rc : cut off radius
X : position of the nucleus
support of the projector function at the center of the box
PAW must not cross the box boundaries !!

fft_paw : performs FFT of f(r)Y_lm(Θ,ϕ) with support B(X, rc)
"""
null(x::Float64,y::Float64,z::Float64) = 0.
null2(x,i,j) = 0.
function fft_paw(rc, X, f, p::params, Npaw; mult = pw_coulomb.multiplier(p.N1), V=null, g=null2)
   @test (X[1]-rc < p.L1) & (X[1]-rc > 0)
   @test (X[2]-rc < p.L2) & (X[2]-rc > 0)
   @test (X[3]-rc < p.L3) & (X[3]-rc > 0)
   Npawtot = paw.Npawtot(Npaw)
   P = zeros(ComplexF64, (2p.N1*mult+1, 2p.N2*mult+1, 2p.N3*mult+1, Npawtot))
   Pout = zeros(ComplexF64, (p.size_psi..., Npawtot))
   FFTW.set_num_threads(4)
   for ipaw in 1:(Npawtot)
      l,n,m = int_to_nl(ipaw,Npaw)
      for i1 in 1:(2*p.N1*mult+1)
         for i2 in 1:(2*p.N2*mult+1)
            for i3 in 1:(2*p.N3*mult+1)
               r = norm(pw_coulomb.coords(p,i1,i2,i3,mult) - X)
               #r = sqrt( ((i1-1)/(2*p.N1*mult+1)*p.L1 - X[1])^2 + ((i2-1)/(2*p.N2*mult+1)*p.L2 -X[2])^2 + ((i3-1)/(2*p.N3*mult+1)*p.L3 - X[3])^2 )
               if (r<rc)
                  P[i1,i2,i3,ipaw] = paw.Y_lm(pw_coulomb.coords(p,i1,i2,i3,mult) - X..., l-1, m)*(f(r,l,n) + g(r,l,n)*V(pw_coulomb.coords(p,i1,i2,i3,mult)...))
               end
            end
         end
      end
      @views P[:,:,:,ipaw] = fft(P[:,:,:,ipaw]) #cannot do on site fft (not doing the right thing)
      @views Pout[:,:,:,ipaw] = pw_coulomb.fft_reshape(P[:,:,:,ipaw], p, mult, 1)
   end
   return Pout
end

function proj_fft(rc, X, p::params, coefpaw::paw.pawcoef; mult = pw_coulomb.multiplier(p.N1))
   f(r,l,n) = polyval(coefpaw.proj[n,l], r)/r
   return fft_paw(rc, X, f, p, coefpaw.Npaw, mult=mult)
end

"""
diff_phi_fft : returns the Fourier coefficients of phi-tilde phi for all the atoms
"""
function diff_phi_fft(rc, X, p::params, coefpaw::paw.pawcoef; mult = pw_coulomb.multiplier(p.N1))
   function f(r,l,n)
      return (paw.R_nl(r, n+l-1, l-1, coefpaw.Z) - paw.tilde_R_nl(r, rc, n+l-1, l-1, coefpaw.Z, coefpaw.tdR[:,n,l]))/r
   end
   return fft_paw(rc, X, f, p, coefpaw.Npaw, mult=mult)
end

"""
Hdiff_phi_fft : returns the Fourier coefficients of H(phi-tilde phi) for all the atoms
"""
function Hdiff_phi_fft(rc, X, p::params, coefpaw::paw.pawcoef, N, Npaw; mult= pw_coulomb.multiplier(p.N1)) #n : atom site
   V(x,y,z) = -coefpaw.Z*sum(1/norm([x,y,z] - X[:,j]) for j in 1:(size(X)[2])) #potential
   Delta_phi(r,n,l) = (paw.E_n(n,coefpaw.Z) + coefpaw.Z/r)*paw.R_nl(r, n, l, coefpaw.Z)/r #(-0.5Δ) ϕ
   tdR(n,l) = Poly(coefpaw.tdR[:,n-l,l+1])(Poly([0,0,1]))
   Delta_tdphi(r,n,l) = paw.C_nl(n,l,coefpaw.Z)*(
   0.5*r^l*polyval(polyder(tdR(n,l),2),r) + (l+1)*r^(l-1.)*polyval(polyder(tdR(n,l)),r)) #0.5Δ \tilde\phi
   function f(r,l,n)
      return Delta_phi(r,n+l-1,l-1) + Delta_tdphi(r,n+l-1,l-1)
   end
   function g(r,l,n)
      return (paw.R_nl(r, n+l-1, l-1, coefpaw.Z) - paw.tilde_R_nl(r, rc, n+l-1, l-1, coefpaw.Z, coefpaw.tdR[:,n,l]))/r
   end
   return fft_paw(rc, X[:,N], f, p, coefpaw.Npaw, mult = mult, V=V, g=g)
end

"""
proj2, diff_phi2, Hdiff_phi2 gives the Fourier coefficients of the PAW functions using direct integration (1D because of the radial symmetry)
int_num works ONLY for l=0 : gives the FFT of f(r)/r*Y_lm(Θ,ϕ)
"""
function int_num(rc, X, f, p::params, Npaw; tol = 1e-10, order = 5)
   Npawtot = paw.Npawtot(Npaw)
   P = zeros(ComplexF64, (p.size_psi..., Npawtot))
   function aux_k(r,k,ipaw,l,n) #K : Fourier mode
      out = r*f(r,l+1,n)*sf_bessel_jl(l,2pi*k*r/p.L1)
      return out
   end
   for ipaw in 1:Npawtot
      l, n, m = int_to_nl(ipaw, Npaw)
      for i3 in 1:(2p.N3+1)
         for i2 in 1:i3
            for i1 in 1:i2
               k_vec = fft_mode.([i1,i2,i3],[p.N1,p.N2,p.N3])
               k = norm(k_vec)
               if k==0.
                  if l>1
                     P[i1,i2,i3,ipaw] = complex(0.)
                  else
                     P[i1,i2,i3,ipaw] = p.Ntot/(p.L1*p.L2*p.L3)*4pi*paw.Y_lm(1.,0.,0.,0,0)*QuadGK.quadgk(r -> r*f(r,l,n), 0, rc, atol=tol, order=order)[1]
                  end
               else
                  aux(r) = aux_k(r,k,ipaw,l-1,n)
                  P[i1,i2,i3,ipaw] = p.Ntot/(p.L1*p.L2*p.L3)*4pi*im^(l-1)*QuadGK.quadgk(aux, 0, rc, atol=tol,rtol=tol,order=order)[1]
                  L = unique(permutations([i1,i2,i3]))
                  for Lj in L[2:end]
                     k_vec = fft_mode.(Lj,[p.N1,p.N2,p.N3])
                     P[Lj...,ipaw] = exp(-2*im*pi*dot(k_vec,X)/p.L1)*paw.Y_lm(-k_vec...,l-1,m)*P[i1,i2,i3,ipaw]
                  end
                  k_vec = fft_mode.([i1,i2,i3],[p.N1,p.N2,p.N3])
                  P[i1,i2,i3,ipaw] = exp(-2*im*pi*dot(k_vec,X)/p.L1)*paw.Y_lm(-k_vec...,l-1,m)*P[i1,i2,i3,ipaw]
               end
            end
         end
      end
   end
   return P
end

function proj_num(rc, X, p::params, coefpaw::paw.pawcoef; tol = 1e-10, order = 5)
   function f(r,l,n)
      if r < rc
         return polyval(coefpaw.proj[n,l], r)
      else
         return 0.
      end
   end
   return int_num(rc, X, f, p, coefpaw.Npaw; tol=tol, order=order)
end

function diff_phi_num(rc, X, p::params, coefpaw::paw.pawcoef; tol = 1e-10, order = 5)
   function f(r,l,n)
      if r < rc
         return paw.R_nl(r, n+l-1, l-1, coefpaw.Z) - paw.tilde_R_nl(r, rc, n+l-1, l-1, coefpaw.Z, coefpaw.tdR[:,n,l])
      else
         return 0.
      end
   end
   return int_num(rc, X, f, p, coefpaw.Npaw; tol=tol, order=order)
end

"""
Hrad_diff_phi_num : computes  the Fourier coefficient of the radial part of H( phi- tilde phi) (direct integration)
"""
function Hrad_diff_phi_num(rc, X, p::params, coefpaw::paw.pawcoef, Npaw; tol = 1e-10, order = 5)
   Hrad_Rnl(r,n,l) = paw.E_n(n,coefpaw.Z)*paw.R_nl(r, n, l, coefpaw.Z) #(-0.5Δ - Z/r) R_nl
   tdR(n,l) = Poly(coefpaw.tdR[:,n-l,l+1])(Poly([0,0,1]))
   Hrad_tdRnl(r,n,l,deriv) = paw.C_nl(n,l,coefpaw.Z)*( 0.5*r^(l+1)*polyval(deriv[2,n-l,l+1],r) + (l+1)*r^l*polyval(deriv[1,n-l,l+1],r)) + coefpaw.Z*paw.tilde_R_nl(r,rc,n,l,coefpaw.Z,coefpaw.tdR[:,n-l,l+1])/r
   #(0.5Δ + Z/r) \tilde R_nl, take into account radial laplacian
   der_tdR = Array{Polynomials.Poly{Float64},3}(undef,2,max(Npaw...),length(Npaw))
   for i in 1:2
      for lpaw in eachindex(Npaw)
         for npaw in 1:Npaw[lpaw]
            der_tdR[i,npaw,lpaw] = polyder(tdR(npaw + lpaw-1, lpaw-1),i)
         end
      end
   end
   function f(r,l,n)
      if r < rc
         return Hrad_Rnl(r,n+l-1,l-1) + Hrad_tdRnl(r,n+l-1,l-1,der_tdR)
      else
         return 0.
      end
   end
#   plot(0:0.01:2rc,[f(x,1,1) for x in 0:0.01:2rc]) : plot OK
   return int_num(rc, X, f, p, coefpaw.Npaw; tol=tol, order=order)
end

"""
Hdiff_phi_pot : computes the Fourier coefficient of the radial part of H(phi-tilde phi) (via a FFT)
"""
function Hdiff_phi_pot(rc, X, p::params, coefpaw::paw.pawcoef, N ; mult=5)
   Y = X[:,1:end .!=N]
   pot(x,y,z) = coefpaw.Z*sum(-1/norm([x,y,z] - Y[:,i3]) for i3 in 1:(length(X[1,:])-1))
   function rad(r,l,n)
      if r < rc
         return paw.C_nl(n,l,coefpaw.Z)*r^l*(GSL.sf_laguerre_n(n-l-1, 2l+1, 2coefpaw.Z*r/n)*exp(-coefpaw.Z*r/n) - polyval(Poly(coefpaw.tdR[:,n-l,l+1]), r^2))
      else
         return 0.
      end
   end
   g(r,l,n) = rad(r,l-1,n+l-1)
   return fft_paw(rc, X[:,N], null2, p, coefpaw.Npaw; mult = mult, V=pot, g=g)
end

function Hdiff_phi_num(rc, X, p::params, coefpaw::paw.pawcoef, n, Npaw ; mult=pw_coulomb.multiplier(p.N1), tol = 1e-14, order = 5) #n : atom site
   diff_phi_mat = Hrad_diff_phi_num(rc, X[:,n], p::params, coefpaw::paw.pawcoef, Npaw; tol=tol, order=order)
   if size(X)[2] == 1
      return diff_phi_mat
   else
      diff_phi_mat += Hdiff_phi_pot(rc, X,  p::params, coefpaw::paw.pawcoef, n ; mult = pw_coulomb.multiplier(p.N1))
      return diff_phi_mat
   end
end

function D_ij(rc, X, coef_PAW, N, Npaw, Z) #N : atom site (working)
   Rnl_diff(r, rc, l1, l2, coef_PAW, j, k) = (paw.R_nl(r, l1+j, l1, Z) - paw.tilde_R_nl(r, rc, l1+j, l1, Z, coef_PAW[:,j,l1+1]))*(paw.R_nl(r, l2+k, l2, Z) - paw.tilde_R_nl(r, rc, l2+k, l2, Z, coef_PAW[:,k,l2+1]))
   DRnl_diff(r, rc, l1, l2, coef_PAW, j, k) = paw.Ddiff(r, rc, l1+j, coef_PAW[:,j,l1+1], l1, Z)*paw.Ddiff(r, rc, l2+k, coef_PAW[:,k,l2+1], l2, Z)
   Npawtot = paw.Npawtot(Npaw)
   D = zeros(ComplexF64,(Npawtot,Npawtot))
   err = zeros(ComplexF64,(Npawtot,Npawtot))
   ind = 0
   for lpaw in eachindex(Npaw)
      l = lpaw -1
      for m in 1:(2lpaw-1)
         for j in 1:Npaw[lpaw]
            for k in 1:Npaw[lpaw]
               kin = QuadGK.quadgk(r -> DRnl_diff(r, rc, l, l, coef_PAW, j, k) + l*(l+1)/r^2*Rnl_diff(r, rc, l, l, coef_PAW, j, k), 0, rc, rtol=1e-10)[1]
               pot1 = QuadGK.quadgk(r -> Z*Rnl_diff(r, rc, l, l, coef_PAW, j, k)/r, 0, rc, rtol=1e-10)[1]
               D[j+ind,k+ind] = 0.5*kin - pot1
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
      function rad(r,n,l)
         if r < rc
            return paw.C_nl(n,l,Z)*r^l*(GSL.sf_laguerre_n(n-l-1, 2l+1, 2Z*r/n)*exp(-Z*r/n) - polyval(Poly(coef_PAW[:,n-l,l+1]), r^2))
         else
            return 0.
         end
      end
      for i in 1:Npawtot
         for j in 1:Npawtot
            l1, n1, m1 = int_to_nl(i, Npaw)
            l2, n2, m2 = int_to_nl(j, Npaw)
            if l1 == l2 == 1 #use radial symmetry to recast to 1D integral
               for i3 in 1:(size(X)[2]-1)
                  R =  norm(X[:,N]-Y[:,i3])
                  pot2, err2 = QuadGK.quadgk(r -> Z/R*Rnl_diff(r, rc, l1-1, l2-1, coef_PAW, n1, n2), 0, rc)
                  D[i,j] += - pot2
                  err[i,j] += err2
               end
            else  #other l1 and l2 values
               realintg(X3) = real( pot(X3...)*rad(norm(X3),n1+l1-1,l1-1)*rad(norm(X3),n2+l2-1,l2-1)*conj(paw.Y_lm(X3...,l1-1,m1))*paw.Y_lm(X3...,l2-1,m2))
               imagintg(X3) = imag( pot(X3...)*rad(norm(X3),n1+l1-1,l1-1)*rad(norm(X3),n2+l2-1,l2-1)*conj(paw.Y_lm(X3...,l1-1,m1))*paw.Y_lm(X3...,l2-1,m2))
               rea, rea_err = int3D(realintg,rc)
               D[i,j] += rea
               err[i,j] += rea_err
               ima, im_err = int3D(imagintg,rc)
               D[i,j] += im*ima
               err[i,j] += im*im_err
            end
         end
      end
      return D, err
   end
end

mutable struct pawfunc
   rc :: Float64 #cut-off radius
   p :: params
   X :: Array{Float64, 2} #positions of the atoms
   Nat :: Integer #number of atoms
   Npaw :: Array{Int64,1} #number of PAW functions
   Npawtot :: Integer
   P #projectors
   Phi #Phi- \tilde\Phi
   Hphi #H(\Phi - \tilde\Phi)
   DH #Array of matrices <phi - \tilde\phi | H  |phi - \tilde\phi >
   DS #Array of matrices <phi - \tilde\phi | phi - \tilde\phi >

   function pawfunc(rc, X::Array{Float64,2}, p::params, Npaw::Array{Int64,1}, Z; proj=proj_fft, diff_phi=diff_phi_fft, Hdiff_phi=Hdiff_phi_fft, args...)
      fpaw = new(rc, p, X, size(X)[2], Npaw, paw.Npawtot(Npaw), nothing)
      Npawtot = paw.Npawtot(Npaw)
      fpaw.DH = zeros(ComplexF64, (Npawtot,Npawtot,fpaw.Nat))
      fpaw.DS = zeros(Float64, (Npawtot,Npawtot,fpaw.Nat))
      fpaw.P = zeros(ComplexF64, (p.size_psi..., Npawtot, fpaw.Nat))
      fpaw.Phi = zeros(ComplexF64, (p.size_psi..., Npawtot, fpaw.Nat))
      fpaw.Hphi = zeros(ComplexF64, (p.size_psi..., Npawtot, fpaw.Nat))
      coefpaw = paw.pawcoef(Z, rc, Npaw; args...)
      for iat in 1:fpaw.Nat
         fpaw.DH[:,:,iat] = D_ij(rc, X, coefpaw.tdR, iat, Npaw, Z)[1]
         fpaw.DS[:,:,iat] = paw.S_ij(rc, Z, coefpaw.tdR, Npaw)
         fpaw.P[:,:,:,:,iat] = proj(rc, X[:,iat], p::params, coefpaw::paw.pawcoef)
         fpaw.Phi[:,:,:,:,iat] = diff_phi(rc, X[:,iat], p::params, coefpaw::paw.pawcoef)
         fpaw.Hphi[:,:,:,:,iat] = Hdiff_phi(rc, X, p::params, coefpaw::paw.pawcoef, iat, Npaw)
      end
      return fpaw
   end
end

"""
tdH_paw computes tilde H tilde psi
"""

function tdH_paw(fpaw::pawfunc, p::params, psi)
   L = p.L1*p.L2*p.L3
   Ppsi = zeros(ComplexF64, (fpaw.Npawtot, fpaw.Nat))
   Hphi_psi = zeros(ComplexF64, (fpaw.Npawtot, fpaw.Nat))
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         Ppsi[ipaw,iat] = 1/p.Ntot^2*L*dot(fpaw.P[:,:,:,ipaw,iat], psi)
         Hphi_psi[ipaw,iat] = 1/p.Ntot^2*L*dot(fpaw.Hphi[:,:,:,ipaw,iat], psi)
      end
   end
   tdpsi = pw_coulomb.ham(p, psi)
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         tdpsi += Ppsi[ipaw,iat]*fpaw.Hphi[:,:,:,ipaw,iat] + Hphi_psi[ipaw,iat]*fpaw.P[:,:,:,ipaw,iat] + (fpaw.DH[:,:,iat]*Ppsi[:,iat])[ipaw]*fpaw.P[:,:,:,ipaw,iat]
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
   Phi_psi = zeros(ComplexF64, (fpaw.Npawtot, fpaw.Nat))
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         Ppsi[ipaw,iat] = 1/p.Ntot^2*L*dot(fpaw.P[:,:,:,ipaw,iat], psi)
         Phi_psi[ipaw,iat] = 1/p.Ntot^2*L*dot(fpaw.Phi[:,:,:,ipaw,iat], psi)
      end
   end
   tdpsi = psi
   for iat in 1:fpaw.Nat
      for ipaw in 1:fpaw.Npawtot
         tdpsi += Ppsi[ipaw,iat]*fpaw.Phi[:,:,:,ipaw,iat] + Phi_psi[ipaw,iat]*fpaw.P[:,:,:,ipaw,iat] + (fpaw.DS[:,:,iat]*Ppsi[:,iat])[ipaw]*fpaw.P[:,:,:,ipaw,iat]
      end
   end
   return tdpsi
end

function energy_vpaw(fpaw::pawfunc, p::params, seed)
   H_1var(psi) = tdH_paw(fpaw, p, reshape(psi, p.size_psi))[:]
   S_1var(psi) = tdS_paw(fpaw, p, reshape(psi, p.size_psi))[:]
   function P(psi)
       meankin = sum(p.kin[i]*abs2(psi[i]) for i = 1:p.Ntot) / (norm(psi)^2)
       return psi ./ (0.1*meankin .+ p.kin[:]) # this should be tuned but works more or less
   end
#   return eigensolvers.eig_lanczos(H_1var, seed[:], B=S_1var, m=5, Imax = 400, do_so=true, norm_A = 6pi^2*p.N1)
   return eigensolvers.eig_pcg(H_1var, seed[:],P=P, B=S_1var, tol=1e-5, maxiter = 400)
end

function tdphi_test(N,L,X,p::params,rc, Z)
   tdphi = zeros(ComplexF64,(2*N+1,2*N+1,2*N+1))
   coef_tdR = paw.coef_tilde_R(rc,1,0,Z)
   for i1 in 1:(2*N+1)
      for i2 in 1:(2*N+1)
         for i3 in 1:(2*N+1)
            r = norm(pw_coulomb.coords(p,i1,i2,i3,1) - X)
            tdphi[i1,i2,i3] = paw.tilde_R_nl(r, rc, 1, 0, Z, coef_tdR)/r
         end
      end
   end
   return fft(tdphi)
end

function test_fft_H(N,L,rc,Npaw,Z)
   X = zeros(3,1)
   X[:,1] = [L/2,L/2,L/2]
   p = params(N,L,X,Z)
   fpaw = pawfunc(rc, X, p, Npaw, Z)
   psi, E, res = energy_vpaw(fpaw, p, tdphi_test(N,L,[L/2,L/2,L/2], p, rc, Z))
   return psi,E
end

function test_num_H(N,L,rc,Npaw,Z; args...)
   X = zeros(3,1)
   X[:,1] = [L/2,L/2,L/2]
   V(r) = -Z/r
   p = params(N,L,X,Z,V)
   fpaw = pawfunc(rc, X, p, Npaw, Z; proj=proj_num, diff_phi=diff_phi_num, Hdiff_phi=Hdiff_phi_num, args...)
   psi, E, res = energy_vpaw(fpaw, p, tdphi_test(N,L,[L/2,L/2,L/2], p, rc,Z))
   return psi, E
end

function guess_H2(rc, X1, X2, p::params, Z)
   tdphi = zeros(ComplexF64,p.size_psi)
   coef_tdR = paw.coef_tilde_R(rc,1,0,Z)
   for i1 in 1:(2*p.N1+1)
      for i2 in 1:(2*p.N2+1)
         for i3 in 1:(2*p.N3+1)
            r1 = norm(pw_coulomb.coords(p,i1,i2,i3,1) - X1)
            r2 = norm(pw_coulomb.coords(p,i1,i2,i3,1) - X2)
            tdphi[i1,i2,i3] = paw.tilde_R_nl(r1, rc, 1, 0, Z, coef_tdR)/r1 + paw.tilde_R_nl(r2, rc, 1, 0, Z, coef_tdR)/r2
         end
      end
   end
   return fft(tdphi)
end

function test_fft_H2(N,L,rc,Npaw,R,Z)
   X1 = [(L-R)/2,L/2,L/2]
   X2 = [(L+R)/2,L/2,L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   V(x,y,z) = -Z/norm([x-(L-R)/2,y-L/2,z-L/2])-Z/norm([x-(L+R)/2,y-L/2,z-L/2])
   p = params(N,L,X,Z)
   fpaw = pawfunc(rc, X, p, Npaw, Z)
   psi, E, res = energy_vpaw(fpaw, p, guess_H2(rc,X1,X2,p,Z))
   return psi, E
end

function test_num_H2(N,L,rc,Npaw,R,Z;args...)
   X1 = [(L-R)/2,L/2,L/2]
   X2 = [(L+R)/2,L/2,L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   V(r) = -Z/r
   p = params(N,L,X,Z,V)
   fpaw = pawfunc(rc, X, p, Npaw, Z; proj=proj_num, diff_phi=diff_phi_num, Hdiff_phi=Hdiff_phi_num, args...)
   psi, E, res = energy_vpaw(fpaw, p, guess_H2(rc,X1,X2,p, Z))
   return psi, E
end

function test_num_H2_cutoff(N,L,rc,Npaw,R,Z;a=0.5*L-R,b=0.5*(L-R),args...)
   @test a<b && L>=3R
   function cut_off(r,a,b)
      if r<=a
         return 1.0
      elseif r<b
         return exp(-(r-a)^6/(b-r)^6)
      else
         return 0.0
      end
   end
   X1 = [(L-R)/2,L/2,L/2]
   X2 = [(L+R)/2,L/2,L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   V(r) = -Z/r*cut_off(r,a,b)
   p = params(N,L,X,Z,V)
   fpaw = pawfunc(rc, X, p, Npaw, Z; proj=proj_num, diff_phi=diff_phi_num, Hdiff_phi=Hdiff_phi_num, args...)
   psi, E, res = energy_vpaw(fpaw, p, guess_H2(rc,X1,X2,p, Z))
   return psi, E
end
end

