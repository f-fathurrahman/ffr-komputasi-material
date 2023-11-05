push!(LOAD_PATH, pwd())

using Polynomials
using Plots
using Cubature
using QuadGK
using GSL
using LinearAlgebra

import VPAWVanderbiltModule
const paw = VPAWVanderbiltModule

# integration on the box [-rc,rc]^3 for functions with singularity at (0,0,0)
# awful compared to radial 1D integration
function int3D(f,rc)
    out = 0
    for i1 in 0:1
       for i2 in 0:1
          for i3 in 0:1
             out += Cubature.hcubature(f, [-rc*(1-i1),-rc*(1-i2),-rc*(1-i3)], [rc*i1,rc*i2,rc*i3],
                abstol=1e-5,reltol=1e-5, maxevals=100000)[1]
          end
       end
    end
    return out
end

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

function TM_test(rc, n, l, Z)
   p = paw.coef_TM(rc, n, l, Z, 1e-8)[1]
   plot(0:0.01:2rc, [paw.V_scr(r, n, l, rc, p, Z) for r in 0:0.01:2rc])
end

function RV_PPtest(rc, n, l, Z)
   p = paw.coef_TM(rc, n, l, Z, 1e-8)[1]
   plot(0:0.01:2rc, [paw.R_V_PP(r, rc, n, l, p, Z) for r in 0:0.01:2rc])
   plot(0:0.01:2rc, [paw.R_nl(r, n, l, Z) for r in 0:0.01:2rc])
end

function td_Rnl_test(rc,n,l,Z)
   p = paw.coef_tilde_R(rc, n, l,Z)
   plot(0:0.001:2rc, [paw.tilde_R_nl(r, rc, n, l, Z, p) for r in 0:0.001:2rc],"b-")
   plot(0:0.001:2rc, [paw.R_nl(r, n, l, Z) for r in 0:0.001:2rc],"r-")
end


function proj_test(rc,n,l,Z,N)
   coef_PAW = zeros(5,N)
   for i in 1:N
      coef_PAW[:,i] = paw.coef_tilde_R(rc,i+l,l,Z)
   end
   coef_PP = paw.coef_TM(rc, l+1, l, Z, 1e-8)[1]
   B = paw.GS_VdB(N, l, Z, rc, coef_PAW, coef_PP)
   plot(0:0.01:2rc, [paw.rad_proj(r, B, n-l, l, Z, rc, coef_PAW, coef_PP, N) for r in 0:0.01:2rc])
end

function proj_test_poly(rc,n,l,Z,N)
   coef_PAW = zeros(5,N)
   for i in 1:N
      coef_PAW[:,i] = paw.coef_tilde_R(rc,i+l,l,Z)
   end
   coef_PP = paw.coef_TM(rc, l+1, l, Z, 1e-8)[1]
   B = paw.GS_VdB(N, l, Z, rc, coef_PAW, coef_PP)
   coef_proj =  paw.coef_rad_proj(B, l, Z, rc, coef_PAW, coef_PP, N)
   plot(0:0.01:rc, [polyval(coef_proj[n-l],r) for r in 0:0.01:rc])
end

function proj_custom_test(rc,n,l,Z,N)
   coef_PAW = zeros(5,N)
   for i in 1:N
      coef_PAW[:,i] = paw.coef_tilde_R(rc,i+l,l,Z)
   end
   coef_PP = paw.coef_TM(rc, l+1, l, Z, 1e-8)[1]
   B = paw.GS_custom(N, l, Z, rc, coef_PAW, coef_PP)
   plot(0:0.01:rc, [paw.rad_proj_custom(r, B, n-l, l, Z, rc, coef_PAW, coef_PP, N) for r in 0:0.01:rc])
end

function proj_custom_test_poly(rc,n,l,Z,N)
   coef_PAW = zeros(5,N)
   for i in 1:N
      coef_PAW[:,i] = paw.coef_tilde_R(rc,i+l,l,Z)
   end
   coef_PP = paw.coef_TM(rc, l+1, l, Z, 1e-8)[1]
   B = paw.GS_custom(N, l, Z, rc, coef_PAW, coef_PP)
   coef_proj =  paw.coef_proj_custom(B, l, Z, rc, coef_PAW, coef_PP, N)
   plot(0:0.01:rc, [polyval(coef_proj[n-l],r) for r in 0:0.01:rc])
end

function ortho_test_rad(rc,Npaw,Z)
   Npawtot = paw.Npawtot(Npaw)
   D = zeros(Float64,(Npawtot,Npawtot))
   coefpaw = paw.pawcoef(Z,rc,Npaw)
   function f(r,n1,l1,m1,n2,l2,m2)
      if r < rc
         return polyval(coefpaw.proj[n1-l1,l1+1],r) * paw.tilde_R_nl(r,rc,n2,l2,Z,coefpaw.tdR[:,n2-l2,l2+1])
      else
         return 0.
      end
   end
   for i1 in 1:Npawtot
      for i2 in 1:Npawtot
         l1,n1,m1 = int_to_nl(i1,Npaw)
         l2,n2,m2 = int_to_nl(i2,Npaw)
         real_f(r) = f(r,n1+l1-1,l1-1,m1,n2+l2-1,l2-1,m2)
         D[i1,i2] = QuadGK.quadgk(real_f, 0.,rc, atol=1e-9, rtol=1e-9)[1]
      end
   end
   return D
end

function ortho_custom_test_rad(rc,Npaw,Z)
   Npawtot = paw.Npawtot(Npaw)
   D = zeros(Float64,(Npawtot,Npawtot))
   coefpaw = paw.pawcoef(Z,rc,Npaw,GS = paw.GS_custom, proj_gen = paw.coef_proj_custom)
   function f(r,n1,l1,m1,n2,l2,m2)
      if r < rc
         return polyval(coefpaw.proj[n1-l1,l1+1],r) * paw.tilde_R_nl(r,rc,n2,l2,Z,coefpaw.tdR[:,n2-l2,l2+1])
      else
         return 0.
      end
   end
   for i1 in 1:Npawtot
      for i2 in 1:Npawtot
         l1,n1,m1 = int_to_nl(i1,Npaw)
         l2,n2,m2 = int_to_nl(i2,Npaw)
         real_f(r) = f(r,n1+l1-1,l1-1,m1,n2+l2-1,l2-1,m2)
         D[i1,i2] = QuadGK.quadgk(real_f, 0.,rc, atol=1e-9, rtol=1e-9)[1]
      end
   end
   return D
end

function ortho_test(rc,Npaw,Z)
   Npawtot = paw.Npawtot(Npaw)
   D = zeros(ComplexF64,(Npawtot,Npawtot))
   coefpaw = paw.pawcoef(Z,rc,Npaw)
   function f(X,n1,l1,m1,n2,l2,m2)
      if norm(X) < rc
         return conj(paw.Y_lm(X...,l1,m1))*polyval(coefpaw.proj[n1-l1,l1+1],norm(X)) * paw.Y_lm(X...,l2,m2)*paw.tilde_R_nl(norm(X),rc,n2,l2,Z,coefpaw.tdR[:,n2-l2,l2+1])/norm(X)^2
      else
         return complex(0.)
      end
   end
   for i1 in 1:Npawtot
      for i2 in 1:Npawtot
         l1,n1,m1 = int_to_nl(i1,Npaw)
         l2,n2,m2 = int_to_nl(i2,Npaw)
         real_f(X) = real(f(X,n1+l1-1,l1-1,m1,n2+l2-1,l2-1,m2))
         imag_f(X) = imag(f(X,n1+l1-1,l1-1,m1,n2+l2-1,l2-1,m2))
         D[i1,i2] = int3D(real_f,rc) + im*int3D(imag_f,rc)
      end
   end
   return D
end

function ortho_test2(rc,Npaw,Z) #other implementation of the projectors \tilde p
   Npawtot = paw.Npawtot(Npaw)
   D = zeros(Complex128,(Npawtot,Npawtot))
   coefpaw = paw.pawcoef(Z,rc,Npaw)
   function f(X,n1,l1,m1,n2,l2,m2,B,coef_PP,N)
      if norm(X) < rc
         return conj(paw.Y_lm(X...,l1,m1))*paw.rad_proj(norm(X), B, n1, l1, rc, coefpaw.tdR[:,:,l1+1], coef_PP, N)* paw.Y_lm(X...,l2,m2)*paw.tilde_R_nl(norm(X),rc,n2,l2,Z,coefpaw.tdR[:,n2-l2,l2+1])/norm(X)^2
      else
         return complex(0.)
      end
   end
   for i1 in 1:Npawtot
      for i2 in 1:Npawtot
         l1,n1,m1 = int_to_nl(i1,Npaw)
         l2,n2,m2 = int_to_nl(i2,Npaw)
         coef_PP = paw.coef_TM(rc, n1+l1-1, l1-1, 1e-8)[1]
         B = paw.GS_VdB(Npaw[l1], l1-1, rc, coefpaw.tdR[:,:,l1], coef_PP)
         N = Npaw[l1]
         real_f(X) = real(f(X,n1+l1-1,l1-1,m1,n2+l2-1,l2-1,m2,B,coef_PP,N))
#         imag_f(X) = imag(f(X,n1+l1-1,l1-1,m1,n2+l2-1,l2-1,m2))
         D[i1,i2] = int3D(real_f,rc)
      end
   end
   return D
end

function pol_proj(rc,n,l)
   L = 0.01:0.01:rc
   coef_PAW = []
   coef_PP = []
   for i in 1:2
      push!(coef_PAW, paw.coef_tilde_R(rc, i, 0))
      push!(coef_PP, paw.coef_TM(rc, i, 0, 1e-8)[1])
   end
   B = paw.GS_VdB(2, 0, rc, coef_PAW, coef_PP)
   isapprox([paw.rad_proj(r, B, n, l, rc, coef_PAW, coef_PP) for r in L], [polyval(paw.coef_rad_proj(B, n, l, rc, coef_PAW, coef_PP)[n], r) for r in L])
end


