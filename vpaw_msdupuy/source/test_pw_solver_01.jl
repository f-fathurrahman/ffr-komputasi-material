using Pkg;
Pkg.activate("../")
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


# Answer for infinite N, L and Z = 1: -0.5,
function H_test(N, L, Z)
    X = zeros(3,1)
    X[:,1] = [L/2, L/2, L/2]
    V(r) = 1/r
    p = PWSolverCoulomb(N, L, X, Z, V)
    psi, E, res = PWCoulombModule.energy(
        p, phi_H(p,[L/2,L/2,L/2],Z), tol=1e-4, maxiter=400
    )
    return psi, E
end

N = 5
L = 10.0
Z = 1.0
H_test(5, 10.0, 1.0)
