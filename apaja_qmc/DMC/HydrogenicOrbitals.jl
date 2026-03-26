
# Hydrogenic orbitals
# ===================
# 1s  (nlm)=(100): ψ1s(r) = (1/√πa0³) * exp(-r/a0)
# 2s  (nlm)=(200): ψ2s(r) = (1/√(32πa0³)) * (2 - r/a0) * eexp(-r/2a0)
# 2px (nlm)=(210) : Ψ2px(r,θ) = (1/√(32πa0³)) * r/a0*exp(-r/2a0)*cosθ = (1/√(32πa0³)) * x/a0 * exp(-r/2a0)
# 2py (nlm)=(211) : Ψ2py(r,θ,ϕ) = (1/√(64πa0³)) * r/a0*exp(-r/2a0)*cosθ*exp(iϕ}  =(real orb.) (1/√(32πa0³)) * y/a0 * exp(-r/2a0)
# 2pz (nlm)=(211) : Ψ2pz(r,θ,ϕ) = (1/√(64πa0³)) * r/a0*exp(-r/2a0)*cosθ*exp(-iϕ} =(real orb.) (1/√(32πa0³)) * z/a0 * exp(-r/2a0)
# a0 = 1 in a.u.

module HydrogenicOrbitals

using GSL  # polynomials
using StaticArrays
import LinearAlgebra: norm
using Printf
using Common

#export phi_1s, phi_2s, phi_2px, phi_2py, phi_2pz, phi_gen, orbital_test
#export psi_1s2s, psi_1s2px, psi_1s2py, psi_1s2pz, psi_2s2pz
export orb, ∇orb, ∇2orb

function init_Beatom_1s2s()
end

function phi_gen(n ::Int64,l ::Int64, m ::Int64, R ::MVector)
    if l>n  
        throw(DomainError())
    end
    if abs(m)>l
        throw(DomainError())
    end
    # normalization
    nor = 1.0
    # nor = sqrt((2/n)^3 * factorial(n-l-1)/(2n*factorial(n+l)))
    # 
    r = norm(R)
    x = R[1]
    y = R[2]
    z = R[3]
    θ = acos(z/r)
    ϕ = atan(y/x)
    ρ = 2r/n
    # these use GSL:
    # radial part
    Rnl = sf_laguerre_n(n-l-1,2*l+1,ρ)*exp(-ρ/2)*ρ^l
    # Spherical harmonic
    # Condon–Shortley phase (-1.0)^m 
    # use abs(m) in Legendre polynomials
    if m==0
        Yml = sf_legendre_Plm(l,0,cos(θ))
    elseif m>0
        Yml = (-1.0)^m*sf_legendre_Plm(l,abs(m),cos(θ))* cos(m*ϕ)
    else
        Yml = (-1.0)^m*sf_legendre_Plm(l,abs(m),cos(θ))* (-1.0)*sin(m*ϕ)
    end
    return nor*Rnl*Yml
end


function phi_2px(x ::MVector)
    r = norm(x)
    return  x[1]*exp(-r/2) 
end

function phi_2py(x ::MVector)
    r = norm(x)
    return  x[2]*exp(-r/2) 
end

function phi_2pz(x ::MVector)
    r = norm(x)
    return  x[3]*exp(-r/2) 
end

function orb(j ::Int64, x ::MVector; Zeff::Float64=1.0) ::Float64
    orb = 0.0 ::Float64
    r = Zeff*norm(x)
    if j==1
        orb = exp(-r)
    elseif j==2
        orb = (2-r)*exp(-r/2) 
    end
    orb
end

function ∇orb(j ::Int64, xvec ::MVector; Zeff::Float64=1.0) :: MVector
    ∇orb = @MVector zeros(3)
    xv = Zeff*xvec
    r = norm(xv)
    
    if j==1
        ∇orb = Zeff*exp(-r)*(-xv/r)
    elseif j==2
        ∇orb = Zeff*exp(-r/2) *(r-4)/2 * xv/r
    else
        println("not done yet")
        exit()
        x = xv[1]
        y = xv[2]
        z = xv[3]
        θ = acos(z./r);
        ϕ = atan(y./x);
    end
    ∇orb
end

function ∇2orb(j ::Int64, xvec ::MVector; Zeff::Float64=1.0) :: Float64
    ∇2orb = 0.0 :: Float64
    xv = Zeff*xvec
    r = norm(xv)
    
    if j==1
        ∇2orb = Zeff^2*exp(-r)*(1-2/r)
    elseif j==2
        ∇2orb = Zeff^2*exp(-r/2)*((6-r)/4 + (r-4)/2*2/r)
    else
        println("not done yet")
        exit()
        x = xv[1]
        y = xv[2]
        z = xv[3]
        θ = acos(z./r);
        ϕ = atan(y./x);
    end
    ∇2orb
end


# some two-electron configurations



@inline function psi_1s2s(x1 ::MVector, x2 ::MVector)
    return phi_1s(x1)*phi_2s(x2)-phi_1s(x2)*phi_2s(x1)
end


@inline function psi_1s2s(x ::MMatrix)
    x1 = x[:,1]
    x2 = x[:,2]
    return phi_1s(x1)*phi_2s(x2)-phi_1s(x2)*phi_2s(x1)
end

@inline function psi_1s2px(x1 ::MVector, x2 ::MVector)
    return phi_1s(x1)*phi_2px(x2)-phi_1s(x2)*phi_2px(x1)
end
@inline function psi_1s2py(x1 ::MVector, x2 ::MVector)
    return phi_1s(x1)*phi_2py(x2)-phi_1s(x2)*phi_2py(x1)
end
@inline function psi_1s2pz(x1 ::MVector, x2 ::MVector)
    return phi_1s(x1)*phi_2pz(x2)-phi_1s(x2)*phi_2pz(x1)
end
    
@inline function psi_2s2pz(x1 ::MVector, x2 ::MVector)
    return phi_2s(x1)*phi_2pz(x2)-phi_2s(x2)*phi_2pz(x1)
end

function orbital_test()
    # assumes orbitals are not normalized
    println("orbital_test:")
    for i in 1:5
        x = @MVector rand(3)
        @printf("phi_1s  %15.10f  %15.10f Δ = %15.10f\n", phi_1s(x), phi_gen(1,0,0,x), phi_1s(x)-phi_gen(1,0,0,x) )
        @printf("phi_2s  %15.10f  %15.10f Δ = %15.10f\n", phi_2s(x), phi_gen(2,0,0,x), phi_2s(x)-phi_gen(2,0,0,x))
        @printf("phi_2px %15.10f  %15.10f Δ = %15.10f\n", phi_2px(x), phi_gen(2,1,1,x), phi_2px(x)-phi_gen(2,1,1,x))
        @printf("phi_2py %15.10f  %15.10f Δ = %15.10f\n", phi_2py(x), phi_gen(2,1,-1,x), phi_2py(x)-phi_gen(2,1,-1,x))
        @printf("phi_2pz %15.10f  %15.10f Δ = %15.10f\n", phi_2pz(x), phi_gen(2,1,0,x), phi_2pz(x)-phi_gen(2,1,0,x))    
    end
    exit()
end
end
