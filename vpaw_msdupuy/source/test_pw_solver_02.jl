push!(LOAD_PATH, pwd())

using FFTW
import PWCoulombModule
using PWCoulombModule: PWSolverCoulomb

function phi_H_real(x, y, z, X, Z)
    dr2 = (x-X[1])^2 +(y-X[2])^2 + (z-X[3])^2
    r = sqrt(dr2)
    return exp(-Z*r)
end

# Solution of the hydrogen atom centered at X
function phi_H(p::PWSolverCoulomb, X, Z)
    mat = zeros(ComplexF64, p.size_psi)
    for i3 in 1:(2*p.N3+1)
        for i2 in 1:(2*p.N2+1)
            for i1 in 1:(2*p.N1+1)
                mat[i1,i2,i3] = phi_H_real(PWCoulombModule.coords(p,i1,i2,i3,1)..., X, Z)
            end
        end
    end
    return fft(mat)
end

function phi_H2(p::PWSolverCoulomb, X1, X2, Z)
    return phi_H(p,X1,Z) + phi_H(p,X2,Z)
end

# X1, X2: coordinates of the atoms
function V_H2(x,y,z, X1, X2)
   r1 = sqrt((x - X1[1])^2 + (y - X1[2])^2 + (z - X1[3])^2)
   r2 = sqrt((x - X2[1])^2 + (y - X2[2])^2 + (z - X2[3])^2)
   return -1/r1 - 1/r2
end

function H2_test(N,L,R,Z)
   X1 = [(L+R)/2, L/2, L/2]
   X2 = [(L-R)/2, L/2, L/2]
   X = zeros(3,2)
   X[:,1] = X1
   X[:,2] = X2
   V(r) = -Z/r
   p = PWSolverCoulomb(N,L,X,Z,V)
   psi, E, res = PWCoulombModule.energy(p, phi_H2(p,X1,X2,Z), tol=1e-10, maxiter=400)
   return psi, E
end

psi, E = H2_test(50, 20, 1.2, 1.0)


