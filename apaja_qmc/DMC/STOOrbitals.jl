module STOOrbitals
using StaticArrays
using ForwardDiff
#import LinearAlgebra: norm

using Common
using GSL  # polynomials

# slightly faster norm 
@inline norm(x) = sqrt(sum(abs2,x))

#=
Here Slater-type atomic orbitals (STO) are
ϕ_j(rvec) = r^{n-1} exp(-ζ_j r) S_{lm}(θ,ϕ) , 


 1s   n=0, l=0, m=0 
 2s   n=1, l=0, m=0
 2p   n=2, l=1, m=-1,0,1
 3s   n=3, l=0, m=0
 3p   n=3, l=1, m=-1,0,1   
in filling order 1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, 4d, 5p, 6s, 4f, 5d, 6p


The first ones are
s:
Y_{00}(θ,ϕ)     = 1/2*sqrt(1/pi)

p:
Y_{1,-1}(θ,ϕ)   = 1/2*sqrt(3/(2pi)) sin(θ) e^{-iϕ}
Y_{10}(θ,ϕ)     = 1/2*sqrt(3/pi) cos(θ) 
Y_{11}(θ,ϕ)     = -1/2*sqrt(3/(2pi)) sin(θ) e^{iϕ}

d:
Y_{2,-2}(θ,ϕ)   = 1/4*sqrt(15/(2pi)) sin^2(θ) e^{-2iϕ}
Y_{2,-1}(θ,ϕ)   = 1/2*sqrt(15/(2pi)) sin(θ) cos(θ) e^{-iϕ}
Y_{20}(θ,ϕ)     = 1/4*sqrt(5/pi) (3 cos^2(θ)-1)
Y_{21}(θ,ϕ)     = -1/2*sqrt(15/(2pi)) sin(θ)cos(θ) e^{iϕ}
Y_{2,-2}(θ,ϕ)    = 1/4*sqrt(15/(2pi)) sin^2(θ) e^{2iϕ}



The real spherical harmonics  S_{lm} (θ,ϕ) are

                  sqrt(2) (-1)^m Re Y_{l|m|}(θ,ϕ) , for m>0
S_{lm}(θ,ϕ) =     Y_{l0}(θ,ϕ)            , for m=0
                  sqrt(2) (-1)^m Im Y_{lm}(θ,ϕ) , for m<0


s:
S_{00}(θ,ϕ)     = 1/2*sqrt(1/pi)                  ϕ_1(rvec) = c_1 exp(-ζ_1 r) 1/2*sqrt(1/pi)

p: (actually px and py are linear combinations of m=1 and m=-1 states, but I arbitrarily assign one m=-1 to px and m=1 to py) 
S_{1,-1}(θ,ϕ)   = 1/2*sqrt(3/pi) sin(θ)sin(ϕ)   =>  2py : r*S_{1,-1}(θ,ϕ) = 1/2*sqrt(3/pi)*y  
S_{10}(θ,ϕ)     = 1/2*sqrt(3/pi) cos(θ)         =>  2pz : r*S_{10}(θ,ϕ)   = 1/2*sqrt(3/pi)*z   
S_{11}(θ,ϕ)     = 1/2*sqrt(3/pi) sin(θ)cos(ϕ)   =>  2px : r*S_{11}(θ,ϕ)   = 1/2*sqrt(3/pi)*x  

d:
S_{2,-2}(θ,ϕ)   = 1/4*sqrt(15/pi) sin^2(θ)sin(2ϕ) 
S_{2,-1}(θ,ϕ)   = 1/4*sqrt(15/pi) sin(2θ)sin(ϕ)
S_{2,0}(θ,ϕ)    = 1/4*sqrt(5/pi) (3cos^2(θ)-1)
S_{2,1}(θ,ϕ)    = 1/4*sqrt(15/pi) sin(2θ)cos(ϕ)
S_{2,2}(θ,ϕ)    = 1/4*sqrt(15/pi) sin(θ)cos(2ϕ)

Atom states are eigenstate of S^2 and S_z, so for p-states etc. take linear combinations,
1s: ϕ_1(rvec) = c_1 * exp(-ζ_1 r) 
2s: ϕ_2(rvec) = c_2 * exp(-ζ_2 r) 
2p: ϕ_3(rvec) = c_3 * exp(-ζ_3 r) r (S_{1,-1}+S_{10}+S_{11})(θ,ϕ)
          = c_3 * exp(-ζ_3 r) * (x+y+z)
ϕ_4(rvec) = c_4 * exp(-ζ_4 r) r^2 (S_{2,-2}+S_{2-1}+S_{20}+S_{21}+S_{22})(θ,ϕ)

=#




# R_nl(r) = N r^{n-1} exp(-ζr)
#  and multiply with S_lm(θ,ϕ)

#The real spherical harmonics S_{lm} (θ,ϕ) are 
#                  sqrt(2) (-1)^m Re Y_{l|m|}(θ,ϕ) , for m>0
#S_{lm}(θ,ϕ) =     Y_{l0}(θ,ϕ)            , for m=0
#                  sqrt(2) (-1)^m Im Y_{lm}(θ,ϕ) , for m<0


#using SphericalHarmonics
#Y = computeYlm(pi/2, pi/3, lmax = 1, m_range = SphericalHarmonics.ZeroTo, SHType = SphericalHarmonics.RealHarmonics())

# Spherical harmonics
# Condon–Shortley phase (-1.0)^m 
# use abs(m) in Legendre polynomials

#=
function Y(l ::Int64, m ::Int64, θ ::T, ϕ ::Z) where T<:Number where  Z<:Number
    if m==0
        Ylm = sf_legendre_Plm(l,0,cos(θ))
    elseif m>0
        Ylm = (-1.0)^m*sf_legendre_Plm(l,abs(m),cos(θ))* exp(im*m*ϕ)
    else
        Ylm = (-1.0)^m*sf_legendre_Plm(l,abs(m),cos(θ))* (-1.0)* exp(-im*m*ϕ)
    end
    return Ylm    
end



@inline function S(l ::Int64, m ::Int64, θ ::T, ϕ ::Z) where T<:Number where  Z<:Number
    if l==0
        S = 1/2*sqrt(1/pi)
    elseif l==1
        if m==-1
            # 2py
            S = 1/2*sqrt(3/pi)*sin(θ)*sin(ϕ)
        elseif m==0
            # 2pz
            S = 1/2*sqrt(3/pi)*cos(θ)
        elseif m==1
            # 2px
            S = 1/2*sqrt(3/pi)*sin(θ)*cos(ϕ)
        else
            error("illegal m")
        end
    elseif l==2
        if m==-2
            S = 1/4*sqrt(15/pi)*sin(θ)^2*sin(2ϕ)
        elseif m==-1
            S =  1/4*sqrt(15/pi)*sin(2θ)*sin(ϕ)
        elseif m==0
            S = 1/4*sqrt(5/pi)*(3cos(θ)^2-1)
        elseif m==1
            S = 1/4*sqrt(15/pi)*sin(2θ)*cos(ϕ)
        elseif m==2
            S = 1/4*sqrt(15/pi)*sin(θ)*cos(2ϕ)
        else
            error("illegal m")
        end
    else
        error("orb l>2 not implemented")
    end
    #=
    if m>0
        S = sqrt(2.0) * (-1)^m * real(Y(l,m,θ,ϕ))
    elseif m==0
        S = Y(l,m,θ,ϕ)
    else
        S = sqrt(2.0) * (-1)^m * imag(Y(l,m,θ,ϕ))
    end
    =#
    return S
end
=#

# normalization of radial STO
@inline norm_STO(n,ζ) = sqrt.( (2*ζ).^(n+1)/factorial(2n) )

# Orbitals
# ========
@inline function orb(rvec ::MVector, ao ::AtomicOrbital)
        
    r = norm(rvec)
    x,y,z = rvec
    n = ao.n
    ζ = ao.exponents
    c = ao.coefficients
    nor = norm_STO(n,ζ)
    
    
    if ao.orbital_type =="1s"
        orb = sum(@. c * nor * exp(-ζ*r))
    elseif ao.orbital_type =="2s"
        orb = sum(@. c * nor * r * exp(-ζ*r))
    elseif ao.orbital_type =="2px"
        orb = sum(@. c * nor * x * exp(-ζ*r))
    elseif ao.orbital_type =="2py"
        orb = sum(@. c * nor * y * exp(-ζ*r))
    elseif ao.orbital_type =="2pz"
        orb = sum(@. c * nor * z * exp(-ζ*r))
    else
        error("orbital not implemented")
    end    
    orb
end




# ∂/∂r R_nl(r) = ((n-1)/r-ζ) R_{nl}(r)
@inline function ∇orb(rvec ::MVector, ao::AtomicOrbital)
    # ∇orb[:,k] is ∇orb kth component
    r = norm(rvec)
    x,y,z = rvec
    n = ao.n
    ζ = ao.exponents
    c = ao.coefficients
    nor = norm_STO(n,ζ)
    cexpo =  @. c * nor * exp(-ζ*r)
    
    if ao.orbital_type =="1s"
        # orb = c * nor * exp(-ζ*r)
        ∇orb = sum(cexpo .* (-ζ) ) * rvec/r
    elseif ao.orbital_type =="2s"
        # orb = c * nor * r * exp(-ζ*r)
        ∇orb = sum(cexpo .* (1 .- ζ.*r)) * rvec/r
    elseif ao.orbital_type =="2px"
        # orb = c * nor * x * exp(-ζ*r)
        ∇orb = sum(cexpo) .* [1,0,0] + sum(cexpo*x .* (-ζ))*rvec/r
    elseif ao.orbital_type =="2py"
        # orb = @. c * nor * y * exp(-ζ*r)
        ∇orb = sum(cexpo) .* [0,1,0] + sum(cexpo*y .* (-ζ))*rvec/r        
    elseif ao.orbital_type =="2pz"
        # orb = @. c * nor * z * exp(-ζ*r)
        ∇orb = sum(cexpo) .* [0,0,1] + sum(cexpo*z .* (-ζ))*rvec/r
    else
        error("orbital not implemented")
    end

    
    #=
    println("-"^80)
    @show ao.orbital_type 
    @show ∇orb
    h = 1.e-6
    nugr = zeros(3)
    f0 = orb(rvec ::MVector, ao ::AtomicOrbital)
    for k in 1:3
        rvec[k] += h
        fp = orb(rvec ::MVector, ao ::AtomicOrbital)
        nugr[k] = (fp-f0)/h
        rvec[k] -= h
    end
    @show nugr
    diff = ∇orb-nugr
    @show diff
    if any(i->i>1e-4, abs.(diff))
        error("deviation")
    end
    =#
    ∇orb
end

# ∇^2 R_{nl}(r) =  (n(n-1)/r^2 - 2nζ/r + ζ^2) R_{nl}(r)
# ∇^2 Y{lm}(θ,ϕ) = -l(l+1)/r^2 Y Y{lm}(θ,ϕ)  
@inline function ∇2orb(rvec ::MVector, ao::AtomicOrbital)
    r = norm(rvec)
    x,y,z = rvec
    θ = acos(z/r);
    ϕ = atan(y/x);

    ζ = ao.exponents
    c = ao.coefficients
    n = ao.n
    l = ao.l
    NSTO =  norm_STO(n,ζ)
    
    if ao.orbital_type =="1s"
        ∇2orb = sum(@. NSTO*c*exp(-ζ*r) *(ζ^2 - ζ*2/r))        
    elseif ao.orbital_type =="2s"
        ∇2orb =  sum(@.  (2/r^2 - 4*ζ/r + ζ^2) * r*NSTO*c*exp(-ζ*r))
    elseif ao.orbital_type =="2px"
        ncexpo = @. NSTO*c*exp(-ζ*r)
        t = x * (ζ.^2 .- 4ζ/r)
        ∇2orb = sum(ncexpo .* t)
    elseif ao.orbital_type =="2py"
        ncexpo = @. NSTO*c*exp(-ζ*r)
        t = y * (ζ.^2 .- 4ζ/r)
        ∇2orb = sum(ncexpo .* t)
    elseif ao.orbital_type =="2pz"
        ncexpo = @. NSTO*c*exp(-ζ*r)
        t = z * (ζ.^2 .- 4ζ/r)
        ∇2orb = sum(ncexpo .* t)
    else
        @show  ao.orbital_type
        error("Not implemented")
    end

    #=
    println("-"^80)
    @show ao.orbital_type 
    @show ∇2orb
    h = 1.e-4
    nu_∇2 = 0.0
    for k in 1:3
        f0 = orb(rvec, ao)
        rvec[k] += h
        fp = orb(rvec, ao)
        rvec[k] -= 2h
        fm = orb(rvec, ao)
        nu_∇2 += (fm- 2*f0 +fp)/(h^2)
        rvec[k] += h
    end
    @show nu_∇2 
    =#
    
    ∇2orb
end



end
