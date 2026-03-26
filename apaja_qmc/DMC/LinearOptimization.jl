# A few linearptimization routines, some copied from QMC_LinearOpt_*.jl

module LinearOptimization
using StaticArrays
using Common
using LinearAlgebra: eigen
using Printf
export solve_Δα!, update_a_opt



H = Matrix{Float64}(undef,1,1)
S = Matrix{Float64}(undef,1,1)
ΨiΨjEL = Matrix{Float64}(undef,1,1)
ΨiELj = Matrix{Float64}(undef,1,1)
ΨiΨj  = Matrix{Float64}(undef,1,1)
Ψi =  Vector{Float64}(undef,1)
ΨiEL = Vector{Float64}(undef,1)
ELi = Vector{Float64}(undef,1)


function get_H_S!(walker, 
                  Ψ_is ::Matrix{Float64},
                  EL_is ::Matrix{Float64},
                  ELs ::Vector{Float64})
    global H, S, ΨiΨjEL, ΨiELj, ΨiΨj, Ψi, ΨiEL, ΨELi, ELi
    npara = size(Ψ_is,1)
    Nw = length(walker)
    
    # averages in Eqs. (53), (54)
    EL = 0.0 
    Ψi .= 0.0
    ΨiEL .= 0.0
    ELi  .= 0.0
    ΨiΨjEL .= 0.0
    ΨiELj  .= 0.0
    ΨiΨj  .= 0.0

            
    for iw ∈ 1:Nw
        # sanity check
        if walker[iw].Ψ == 0.0
            error("LinearOptimization walker[iw].Ψ=0; set the walker trial wave function value")
        end

        Ψ_i_R = Ψ_is[:,iw]
        EL_i_R = EL_is[:,iw]
        EL_R = ELs[iw]
        Ψ_0_R = walker[iw].Ψ

        # add to averages:
        @inbounds for i = 1:npara
            Ψi[i] +=  Ψ_i_R[i]/Ψ_0_R           # <Ψ_i/Ψ_0>
            ΨiEL[i] +=  Ψ_i_R[i]/Ψ_0_R * EL_R  # <Ψ_i/Ψ_0 EL>
            ELi[i] += EL_i_R[i]                # <EL_i>
            @inbounds for j = 1:npara
                ΨiΨjEL[i,j] += Ψ_i_R[i]*Ψ_i_R[j]/Ψ_0_R^2 * EL_R  # <Ψ_i/Ψ_0 Ψ_j/Ψ_0 EL>
                ΨiELj[i,j] +=  Ψ_i_R[i]/Ψ_0_R * EL_i_R[j]        # <Ψ_i/Ψ_0 EL_j>
                ΨiΨj[i,j] += Ψ_i_R[i]*Ψ_i_R[j] /Ψ_0_R^2          # <Ψ_i/Ψ_0 Ψ_j/Ψ_0>, may be used in N_i 
            end
        end
        #
        EL += EL_R
        
    end
    mul = 1/Nw
    
    ΨiΨjEL *= mul
    ΨiΨj *= mul
    ΨiELj *= mul        
    Ψi *= mul
    ΨiEL *= mul
    ELi *= mul
    EL *= mul
    # Toulouse-Umrigar J. Chem. Phys. 2008, Eqs. (22), (23), and (24)
    H .= 0.0
    S .= 0.0
    H[1,1] = EL
    S[1,1] = 1.0 # (53b)
    @inbounds for i = 1:npara
        H[1,i+1] = ΨiEL[i] - Ψi[i]*EL + ELi[i] # g_R/2 in (22) or (54c) ; this is the first *row*
        H[i+1,1] = ΨiEL[i] - Ψi[i]*EL # g_L/2 in (22), or (54b) ; this is the first *column*

        @inbounds for j = 1:npara
            # Toulouse-Umrigar J. Chem. Phys. 2008, Eqs. (25), and (26) ; or (54d) and (53c)
            H[i+1,j+1] = ΨiΨjEL[i,j] - Ψi[i]*ΨiEL[j] - Ψi[j]*ΨiEL[i] +
                Ψi[i]*Ψi[j]*EL + ΨiELj[i,j] - Ψi[i]*ELi[j]
            S[i+1,j+1] = ΨiΨj[i,j] - Ψi[i]*Ψi[j]
        end
    end
end



function solve_Δα!(walker, Ψ_is ::Matrix{Float64}, EL_is ::Matrix{Float64}, ELs ::Vector{Float64},
                   a_opt ::Float64, npara ::Int64, nonlin ::Vector{Int64}, Δα:: MVector)

    global H, S, ΨiΨjEL, ΨiELj, ΨiΨj, Ψi, ΨiEL, ΨELi, ELi
    if size(H,1)!=npara+1

        H = Matrix{Float64}(undef,npara+1,npara+1)
        S = Matrix{Float64}(undef,npara+1,npara+1)
        ΨiΨjEL = Matrix{Float64}(undef,npara,npara)
        ΨiELj = Matrix{Float64}(undef,npara,npara)
        ΨiΨj  = Matrix{Float64}(undef,npara,npara)
        Ψi =  Vector{Float64}(undef,npara)
        ΨiEL = Vector{Float64}(undef,npara)
        ELi = Vector{Float64}(undef,npara)
    end

    get_H_S!(walker, Ψ_is, EL_is, ELs)
    
    # stabilize (could also stabilize nonlinear params more)
    for i ∈ 2:npara+1
        H[i,i] += a_opt
    end
    # same as H[diagind(H)[2:npara+1] += a_opt
    #
    #
    # This is sometimes left out, but it stabilized the generalized eigenvalue problem for large a_opt 
    
    #if a_opt>1000.0
    #    for i ∈ 2:npara+1
    #        S[i,i] += a_opt
    #    end
    #end
    
    
    F = eigen(H, S) 
    eigval = real(F.values)
    
    k = argmin(abs.(eigval .- H[1,1])) # "physically reasonable eigenvalue"
    # eigen: "The kth eigenvector can be obtained from the slice F.vectors[:, k]."
    eigv = real(F.vectors)
    #     
    # parameter changes, keep also the 1st element
    Δα_all = eigv[:,k] # this is c*[1,Δα], but scale Δα to get better convergence for nonlinear parameters
    
    # Toulouse-Umrigar 2007
    # ======================
    sc = Δα_all[1]
    # linear parameters
    for i ∈ 1:npara
        # skip nonlinear parameters
        if i ∈ nonlin
            continue
        end
        sc -=  Δα_all[i+1] * Ψi[i]
    end
    # loops do the same as
    #  sc -= sum(Δα_all[i+1] * Ψi[i] for i in 1:npara if i ∉ nonlin)
    # or
    # lin = [i for i ∈ 1:npara if i ∉ nonlin]
    #  sc -= sum(Δα_all[lin.+1] .* Ψi[lin])
    # 
    # so far sc = Δα_all[1] - ∑_i^{linear}  Δα_all[i+1] * Ψi[i]
    #
    # nonlinear parameters
    # 
    # need only one sum 
    Ψ_lin_norm = 0.0
    for j in nonlin        
        for k in nonlin
            # note: shift indices by one
            Ψ_lin_norm += S[j+1,k+1] * Δα_all[j+1] * Δα_all[k+1]
        end
    end
    Ψ_lin_norm /=  Δα_all[1]^2 # normalize Δα_all first component to 1
    #
    # ξ=1 means orthogonalization to Ψ_0
    # ξ=0 means orthogonalization to Ψ_lin
    # general:  orthogonalization to ξ*Ψ_0 + (1-ξ)*Ψ_lin/||Ψ_lin||
    ξ =  0.5 
    sc += (1.0-ξ)*Ψ_lin_norm/( (1.0-ξ) + ξ*sqrt(1.0 + Ψ_lin_norm) )
    
    # normalize solution to [1,Δα], pick Δα    
    Δα .= Δα_all[2:end]/sc
end
end

