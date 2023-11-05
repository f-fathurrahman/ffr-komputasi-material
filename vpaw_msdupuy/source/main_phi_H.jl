using Test
using FFTW
using LinearAlgebra
using FastConv
using DelimitedFiles


# Coulomb potential with cut-off to get smoothness at boundary of the box
function H2_test_cutoff(N,L,R,Z;a=0.5*L-R,b=0.5*(L-R))
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
   X1 = [(L+R)/2,L/2,L/2]
   X2 = [(L-R)/2,L/2,L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   V(r) = -Z/r*cut_off(r,a,b)
   p = pw_coulomb.params(N,L,X,Z,V)
   psi, E, res = pw_coulomb.energy(p,phi_H2(p,X1,X2,Z), tol=1e-10, maxiter=400)
   return psi, E
end

function conv_H2(X1, X2, Nmin, Nmax, L)
   Ns = Nmin:5:Nmax
   Es = zeros(length(Ns))
   Psi = Array(Any, length(Ns))
   for (i,N) in enumerate(Ns)
      println(n)
      psi, E, res = H2_test(N,L,X1,X2)
      Es[i] = E
      # Psi[Integer(round((n-Nmin)/5))+1] = psi
   end
   writedlm("Direct_Energy_45.txt", Es)
   # writedlm("Direct_Psi_45.txt", Psi)
end

function herm(N,L,Z)
   V(x,y,z) = 0.
   p = pw_coulomb.params(N,L,V)
   psi = phi_H(p, [L/2,L/2,L/2], Z)
   pot = zeros(ComplexF64, (4p.N1+1,4p.N2+1,4p.N3+1))
   for i1 in 4p.N1+1
      for i2 in 4p.N2+1
         for i3 in 4p.N3+1
            k_vec = pw_coulomb.fft_mode.([i1,i2,i3],[p.N1,p.N2,p.N3])
            if k_vec == [0,0,0]
               pot[i1,i2,i3] = complex(0.)
            else
               pot[i1,i2,i3] += - p.L1^2*exp(-2*im*pi*dot(k_vec,[L/2,L/2,L/2])/p.L1)/pi/dot(k_vec,k_vec)
            end
         end
      end
   end
   pot_psi = convn(fftshift(pot),fftshift(psi))
   return dot(psi,ifftshift(pot_psi[2p.N1+1:4p.N1+1,2p.N2+1:4p.N2+1,2p.N3+1:4p.N3+1]))
end


function test_V_symmetric()
    X = zeros(3,1)
    X[:,1] = [2.5,2.5,2.5]
    p = params(5,5.,X,3.)
    p.V_grid = fftshift(p.V_grid)
    V_transpose = zeros(ComplexF64,(4p.N1+1,4p.N1+1,4p.N1+1))
    for i1 in 1:4p.N1+1
        for i2 in 1:4p.N2+1
            for i3 in 1:4p.N3+1
                V_transpose[i1,i2,i3] = p.V_grid[4p.N1+2-i1,4p.N2+2-i2,4p.N3+2-i3]
            end
        end
    end
    return norm(p.V_grid - conj.(V_transpose),Inf)
end

# #test for V = -Z/|x| and check if the solution is radial
# function energy_check(V, N1, N2, N3, L1, L2, L3, seed)
#   psi, E, res = energy(V, N1, N2, N3, L1, L2, L3, seed)
#   psi = reshape(psi, (2*N1+1, 2*N2+1, 2*N3+1))
#   ifft!(psi)
#   psi = reshape(psi, (2*N1+1)*(2*N2+1)*(2*N3+1))
#   x = Array(Any, (2*N1+1,2*N2+1,2*N3+1))
#   for i3 in 1:(2*N1+1)
#     for i2 in 1:(2*N2+1)
#       for i1 in 1:(2*N3+1)
#         x[i1,i2,i3] = sqrt(((i1-1)/(2*N1+1)*L1 - L1/2)^2 +
#          ((i2-1)/(2*N2+1)*L2 - L2/2)^2 + ((i3-1)/(2*N3+1)*L3 - L3/2)^2  )
#       end
#     end
#   end
#   return reshape(x, (2*N1+1)*(2*N2+1)*(2*N3+1)), psi
# end


