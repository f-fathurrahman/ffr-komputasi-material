
using Printf
import Distributions: Normal
using StaticArrays
using Statistics

import LinearAlgebra: ⋅
# slightly faster norm 
@inline norm(x) = sqrt(sum(abs2,x))


# local modules:
push!(LOAD_PATH,".")
using Common
using Utilities
using QMC_Statistics
using SlaterDeterminants
using STOOrbitals
using Jastrow
using VMCstep

# for testing purposes only:
#using Random
#Random.seed!(123456)



possible_args=["atom","order","tau","basis_set"]
arg_dict = argparse(possible_args)
# use command line values or defaults
atom = get(arg_dict, "atom", "")
order_in = get(arg_dict, "order", 1) 
τ_in = get(arg_dict, "tau", 0.001)
basis_set = get(arg_dict, "basis_set", "STO")


const order = order_in :: Int64
const τ = τ_in :: Float64

const κ = 10.0 :: Float64 # may need adjustment
const accuracy_goal = 0.0001

const Nwx = 20000 :: Int64
const Nw_target = 500 :: Int64

const blocksize = floor(Int64,1/(5*τ)) :: Int64 # just an approx



atom_data = parse_atom_data(atom)

# Print atom data
println("Atom: ", atom_data.atom)
println("Number of electrons (N): ", atom_data.N)
println("Number of spin-up electrons (Nup): ", atom_data.Nup)
println("Number of spin-down electrons (Ndo): ", atom_data.Ndo)
println("Atomic number (Z): ", atom_data.Z)
println("Exact energy (Eexact): ", atom_data.Eexact)

const N = atom_data.N
const Z = atom_data.Z
const Nup =atom_data.Nup
const Ndo = atom_data.Ndo
const Eexact = atom_data.Eexact



# cutoff to prevent uncontrolled branching near nodal boundaries
# ZSGMA by Zen, Sorella, Gillan, Michaelides and Alfe, PRB 93, 241118 (2016) Eq. (5) and (6)
# factor 0.2 is just a pick 
const Ecut = 0.2*sqrt(N/τ) ::Float64

jpara = parse_Jastrow_data(atom, basis_set)

println("Jastrow parameters:")
@show jpara.param


# Select atomic orbital functions
SlaterDeterminants.init(basis_set)
all_aos = parse_basis_set_data(atom, basis_set)

# Read lc coefficients from file  
cnls = parse_coeffients_data(atom, basis_set)

# create linear combination of ao's 
lc, npara_, nonlin_, lc_ups_, lc_dos_ = set_linear_combination_of_orbitals(atom_data, all_aos, cnls)


# fix lists of orbitals
const lc_ups = lc_ups_
const lc_dos = lc_dos_

# sanity check
for aos in lc.orbitals
    if length(aos) != N
        @show length(aos), N
        error("lc.orbitals must have length N to accomdate all elecrons")
    end
end




suffix_ = string(PROGRAM_FILE,"_",atom,"_order=",order,"_tau=",τ,"_basis_set=",basis_set)

const suffix = suffix_ ::String
const Efile = string("E.",suffix)
const Eallfile = string("E.tmp.",suffix)
rm(Eallfile,force=true)


println("="^40)
@show suffix 
println("Energy output file: ",Efile)
println("="^40)


# collect parameters to one struct
const params = WfParams(lc, jpara)


SlaterDeterminants.init(basis_set)
# 2 for parallel spins, 4 for antiparallel
const s = spinfactor(Nup, Ndo)



function generate_walkers(Nw ::Int64, Nwx ::Int64, ϕ ::Function)
    println("generating $Nw walkers, allocating $Nwx walkers")
    r = [MMatrix{dim,N}(rand(dim,N) .- 0.5) for i in 1:Nwx]

    
    pdata = nothing
    walker = [Walker(;R=r[i], Ψ = ϕ(r[i]), E=0.0, alive=0, age=0.0, weight=1.0, pure_data=pdata) for i in 1:Nwx]

    vmc_params = VMC_Params(0, 0, 1.0)
    
    println("thermalizing walkers")
    E = 0.0
    alives = Set(1:Nw)  # set of alive walkers
    for iw in 1:Nw
        for i in 1:100
            vmc_step!(walker[iw].R, vmc_params, x -> ϕ(x))                
            if i%10 == 0
                adjust_step!(vmc_params)
            end
        end
        walker[iw].Ψ = ϕ(walker[iw].R)
        walker[iw].E = EL(walker[iw].R, ϕ, params)
        walker[iw].alive = 1
        E += walker[iw].E

    end
        
    for p in vmc_params
        @show p.step
    end
    E /= Nw
    println("<E> = $E")
    return walker, E, alives
end

# version with default parameters
@inline @fastmath function φ_T(R ::MMatrix{dim,N,Float64}) 
    return φ_T(R, params)
end



@inline function φ_T(R ::MMatrix{dim,N,Float64}, params ::WfParams)
    return φ_T(R, params.lc, params.jpara)
end

@inline function φ_T(R ::MMatrix{dim,N,Float64}, lc ::LinearCombination, jpara ::JastrowParameters) 
    DD = 0.0    
    for (i,aos) in enumerate(lc.orbitals)        
        ups = aos[lc_ups[i]]
        dos = aos[lc_dos[i]]
        DD += lc.coeffs[i] * D_up(R, ups) * D_do(R, dos)
    end
    if isapprox(DD,0.0)
        error("φ_T: DD = 0, probably bad orbitals")
    end
    φ_T = φ_T_J(R, jpara, s) * DD
end


@inline function EL(R ::MMatrix{dim,N,Float64}, params::WfParams)
    return EL(R, params.lc, params.jpara)
end

@inline function EL(R ::MMatrix{dim,N,Float64}, lc ::LinearCombination, jpara ::JastrowParameters)
    
    VL = 0.0 

    ∇DD = zeros(dim,N)
    TL = 0.0
    DDsum = 0.0
    for (i,aos) in enumerate(lc.orbitals)
        ups = aos[lc_ups[i]]
        dos = aos[lc_dos[i]]
        # see comments before F(R)
        cDD = lc.coeffs[i] * D_up(R, ups) * D_do(R, dos)
        TL += -D* cDD * ( ∇2_up(R, ups) + ∇2_do(R, dos) )        
        ∇DD +=  cDD *( F_up(R, ups) + F_do(R, dos) )
        DDsum +=  cDD
    end
    TL /= DDsum
    ∇DD /= (DDsum*2)    
    ∇J   = F_J(R, jpara, s)/2     # (∇ e^J)/e^J = ∇J = F/2
    
    TL += -D*2* ∇DD ⋅ ∇J
    
    α, β1, β2, α12 = jpara.param
    ∇2J  = 0.0
    for k in 1:N
        rk  = norm(R[:,k]) 
        ∇2J += -α*2/rk #  2 is dim-1
        VL += -Z/rk
        for i in k+1:N
            ski = s[k,i]
            rki_vec = R[:,k]-R[:,i] 
            rki = norm(rki_vec) 
            VL += 1/rki
            if ski==2
                ∇2J += 4* α12/(ski*(1+β1*rki)^3)/rki
            else
                ∇2J += 4* α12/(ski*(1+β2*rki)^3)/rki
            end
        end        
    end    
    TL += -D * ( sum(∇J.^2) + ∇2J )
    EL = TL + VL
    EL    
end


# Total drift  F = 2*∇φ_T/φ_T
# F_up returns 2∇D_up/D_up for a single determinant, same with F_do
# so for a linear combination of two determinants the Slater part is 
# F =  2*∇(D_up(1)D_do(1) + D_up(2)D_do(2))/ (D_up(1)D_do(1) + D_up(2)D_do(2))
# = (2*∇D_up(1)*D_do(1) + 2*∇D_do(1)*D_up(1)+ 2*∇D_up(2)*D_do(2) + 2*∇D_do(2)*D_up(2)) /(D_up(1)D_do(1) + D_up(2)D_do(2))
# = (F_up(1)*D_up(1)*D_do(1) + F_do(1)*D_do(1)*D_up(1)+ F_up(2)*D_up(2)*D_do(2) + F_do(2)*D_do(2)*D_up(2))
#        /(D_up(1)D_do(1) + D_up(2)D_do(2))
# If there's only one determinant, this simplifies to
# (F_up(1)*D_up(1)*D_do(1) + F_do(1)*D_do(1)*D_up(1))/(D_up(1)D_do(1)) = F_up(1) + F_do(1).

@inline function get_F!(R ::MMatrix{dim,N,Float64}, F ::MMatrix{dim,N,Float64} ; UNR=true)
    F .= 0.0
    DDsum = 0.0
    for (i,aos) in enumerate(lc.orbitals)
        ups = aos[lc_ups[i]]
        dos = aos[lc_dos[i]]
        #
        cDD = lc.coeffs[i] * D_up(R, ups) * D_do(R, dos)
        F .+= cDD * ( F_up(R, ups) + F_do(R, dos) )
        DDsum += cDD
    end
    F ./= DDsum  # blows up if DDsum≈0.0
    F .+= F_J(R, jpara, s)
    #
    if UNR
        UNR_scaling!(F)
    end
    #
end

# wrapper for drift testing *without* UNR_scaling
function F_wrapper(R ::MMatrix{dim,N,Float64})
    F = @MMatrix zeros(dim,N)
    get_F!(R, F; UNR=false)
    return F
end


# Drift near nodal surfaces gets very large, which can be approximately
# corrected as in (UNR) Umrigar, Nightingale, Runge, JChemPhys 99, 2865 (1993).
@inline function UNR_scaling!(F ::MMatrix{dim,N,Float64})
    for i = 1:N
        vsq = sum(F[:,i].^2)/4   # v = F/2
        sc = (-1.0 + sqrt(1.0 + 2.0*vsq*τ))/(vsq*τ) 
        # vsq may be very small, no need to scale
        # sc -> 1 as τ -> 0
        sc = min(1.0, sc) 
        F[:,i] = sc*F[:,i]
    end
end



@inline function lnG(xp ::MMatrix{dim,N,Float64}, x ::MMatrix{dim,N,Float64}, F ::MMatrix{dim,N,Float64})         
    # only parts not symmetric in x'<->x
    # G(x'<-x,τ) = exp(-(x'-x-dτF(x))^2/(4Dτ)
    #     
    get_F!(x, F)
    # re-use F
    @. F = xp - x - D*τ*F
    lnG = -sum(F.^2)/(4*D*τ) 
end


# walker update
@inline function diffusion_drift_step!(x ::MMatrix{dim,N,Float64},
                                       F ::MMatrix{dim,N,Float64},
                                       y ::MMatrix{dim,N,Float64},
                                       z ::MMatrix{dim,N,Float64})
    #
    # help arrays y and z
    #
    if order==1
        # diffusion(τ)+drift(τ)
        η = reshape(rand(Normal(0,1),dim*N),(dim,N))
        get_F!(x, F)
        @. x += sqrt(2*D*τ)*η + D*τ*F 
        
    elseif order==2
        # dangerous! diffusion may have taken near nodal surface
        # UNR scaling can't save if F = ***/0.0  
        
        # diffusion(τ/2)+drift(τ)+diffusion(τ/2)
        #
        # step 1)
        η = reshape(rand(Normal(0,1),dim*N),(dim,N))        
        @. y = x + sqrt(D*τ)*η
        get_F!(y, F, UNR=true)  
        # step 2)
        @. z = y + D*τ/2*F
        get_F!(z, F, UNR=true)
        # step 3)
        η = reshape(rand(Normal(0,1),dim*N),(dim,N))
        @. x = y + D*τ*F + sqrt(D*τ)*η
        
        #
        #= *** check validity ***
        # drift(τ/2)+diffusion(τ)+drift(τ/2)
        #
        # Drift(τ/2) 2nd order accurately:
        get_F!(x, F, UNR=true)
        @. y = x + D*τ/4*F
        get_F!(y, F, UNR=true)
        @. z = x + D*τ/2*F
        # diffusion(τ):
        η = reshape(rand(Normal(0,1),dim*N),(dim,N))
        @. x = z + sqrt(2*D*τ)*η
        # Drift(τ/2) 2nd order accurately:
        get_F!(x, F, UNR=true)
        @. y = x + D*τ/4*F
        get_F!(y, F, UNR=true)
        @. x += D*τ/2*F
        =#
    elseif order==3
        dt = 2*D*τ
        eta = reshape(rand(Normal(0,1),dim*N),(dim,N))
        @. x = x + sqrt(dt/2)*eta
        get_F!(x, F)
        @. y = x + dt/2*F/2
        get_F!(y, F)
        eta = reshape(rand(Normal(0,1),dim*N),(dim,N)) 
        @. x += dt*F2/2 + sqrt(dt/2)*eta
    end
end

function run_VMC!(walker ::Vector{Walker}, vmc_params ::VMC_Params, Eallfile ::String) ::Float64   
    Nw = length(walker)
    println("thermalizing $Nw walkers using VMC")
    @inbounds for iw in 1:Nw
        R = walker[iw].R
        @inbounds for i in 1:50
            vmc_step!(R, vmc_params, φ_T)
            if i%10 == 0
                adjust_step!(vmc_params)
            end
        end
    end
    println("VMC       step = $(vmc_params.step)")
    println("VMC acceptance = $(vmc_params.Naccept*100.0/vmc_params.Ntry)")

    #
    # VMC
    #
    Estat = init_stat(1, 10)
    E_ave = 0.0
    while true
        E_ave *= 0.0 # make sure we use the same E_ave as outside 
        @inbounds for iw in 1:Nw
            R = walker[iw].R
            vmc_step!(R, vmc_params, φ_T)
            E = EL(R, params)
            walker[iw].E = E    
            walker[iw].Ψ = φ_T(R)
            E_ave += E
        end
        open(Eallfile,"a") do f
            println(f,E_ave/Nw)
        end
        add_sample!(Estat, E_ave/Nw)
        adjust_step!(vmc_params)
        if Estat.finished            
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)
            @printf("VMC <E> = %.10f +/- %.10f \n", Eb_ave, Eb_std)
            # short VMC 
            if (Eb_std<0.001 && Nb>5) || Nb==20
                output_MCresult(Eb_ave, Eb_std)
                return Eb_ave 
            end
            # just VMC 
            #if (Eb_std<0.0001 && Nb>5) 
            #    output_MCresult(Eb_ave, Eb_std)
            #    output_MCresult(Eb_ave, Eb_std)
            #    println("VMC done")
            #    exit()
            #end
        end        
    end
end

function main()
    #
    # Main program 
    #
   

    # ------------------
    # preallocate arrays
    # drift
    F = @MMatrix zeros(dim,N)
    # help arrays
    help1 = @MMatrix zeros(dim,N)
    help2 = @MMatrix zeros(dim,N)
    
    #
    # -----------------
    # Checks
    
    Nw = 10
    walker =  init_walkers(Nw, Nw, dim, N, φ_T=φ_T, scale=1.0)
    
    
    for i in 1:Nw
        num_check_∇S(walker[i].R, x->φ_T_J(x, jpara, s), x->F_J(x, jpara, s))
    end
    println("Jastrow drift check passed")


    for i in 1:Nw
        num_check_∇S(walker[i].R, φ_T, F_wrapper)
    end
    println("Full drift check passed")

    println("Checking numerically EL against φ_T")

    
    for i in 1:Nw
        R = walker[i].R
        E = EL(R, params)
        num_check_EL(R, E, φ_T, x->Common.V(x,N,Z))
    end
    println("Full EL check passed")
    
    println("checks passed")
    # ---------------------
    #
    
        
    # init walkers
    Nw = Nw_target
    walker = init_walkers(Nw, Nwx, dim, N, φ_T=φ_T, scale=1.0)
    #
    # VMC
    #
    vmc_params = VMC_Params(0, 0, 3.)
    @time E_ave = run_VMC!(walker[1:Nw], vmc_params, Eallfile)
    #exit() # stop after VMC
    #
    # DMC
    #
    ET = E_ave # trial energy, to be updated
    println("ET = $ET")

    alives = Set(1:Nw)

    Ntherm = floor(Int64,2.0/τ) 
    idmc = -Ntherm
    nE = 1
    
    copies = @MVector zeros(Int64,Nwx)
    Ntry = 0 ::Int64
    Nacc = 0 ::Int64
   
    
    # init E measurement
    Estat = init_stat(1, blocksize)
    ETstat = init_stat(1, blocksize)

    while true
        #
        # take a DMC step in each walker
        #
        
        @inbounds for iw in alives
            ELold = walker[iw].E
            Rnew = walker[iw].R # alias
            Rold  = copy(Rnew) # copy
            φTold = walker[iw].Ψ
            #
            diffusion_drift_step!(Rnew, F, help1, help2)
            #
            Ntry +=1
            accept = true
            
            # check we are still in the same nodal cell

            φTnew = φ_T(Rnew)            
            if φTnew < 0
                # reject the move
                accept = false
                # kill walker and skip the rest (not recommended)
                #copies[iw] = 0
                #continue
            end
           
           
            if accept             
                if order==1
                    # Metropolis-Hastings Detailed Balance          
                    Wnew = lnG(Rold, Rnew, F) + 2*log(φTnew) # lnG(new->old) + ln(ϕ_T^2(new))
                    Wold = lnG(Rnew, Rold, F) + 2*log(φTold) # lnG(old->new) + ln(ϕ_T^2(old))
                    if walker[iw].age > 50.0
                        Wnew += (walker[iw].age-50.0)*0.095 # avoid persistent walkers
                    end
                    accept = metro(Wold,Wnew) # takes log weights
                end
            end           
            
            
            if accept
                Nacc +=1
                ELnew = EL(Rnew, params)    
                walker[iw].E = ELnew
                walker[iw].Ψ = φTnew
                walker[iw].age = 0.0
            else
                ELnew = ELold
                Rnew .= Rold
                walker[iw].age += 1.0
            end
            Ebest = E_ave/nE
            Emid = (ELnew+ELold)/2
            # Remedy for singular EL near nodal surface
            Sbar = ET - (Ebest + sign(Emid-Ebest)*min(Ecut, abs(Emid-Ebest)))
            weight = exp(τ*Sbar)
            # basic branching weights that become poor near nodal surface:
            #weight  = exp(-τ*(0.5*(ELold+ELnew)-ET)) ::Float64
            #weight  = exp(-τ*(ELold-ET)) ::Float64

            copies[iw]  = floor(Int64,weight + rand())
        end
        
        # Branching        
        branch!(walker, copies, alives)
        
        # collect energy from walkers
        # note: do this *after* drift, diffusion and branching ! 
        E   = 0.0  
        n   = 0
        max_age = 0.0
        @inbounds for iw in 1:Nwx
            if walker[iw].alive == 1
                max_age = max(max_age, walker[iw].age)
                E += walker[iw].E
                n += 1                
            end
        end
        if max_age>50.0
            @show max_age
        end
        # n is Nw
        E_ave += E/n
        nE += 1
        
        idmc +=1
        
        # update trial energy
        ET = E_ave/nE + κ*log(Nw_target/Nw)

        #E_ave += ET
        #nE += 1
        
        # add new energy data
        if idmc>=0
            add_sample!(Estat, E/n)
            add_sample!(ETstat, ET)
        end
        

        if idmc<0
            @printf("DMC %10d  τ = %.5e <E> = %.10f E = %.10f ET = %.10f <E>_exact = %.10f  %6d Walkers  acceptance = %.5f \n",
                    idmc, τ, E_ave/nE, E/n, ET, Eexact, Nw, 100.0*Nacc/Ntry)
        end
        
        open(Eallfile,"a") do f
            println(f,E/n," ",E_ave/nE," " ,ET," ",Nw)
        end
        # block data
        # ==========
        # screen and file output
        if Estat.finished            
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)
            ETb_ave, ETb_std, E_inputvar2, Nb_ = get_stats(ETstat)

            @printf("DMC %10d τ = %.5e  <E> = %.10f +/- %.10f <ET> = %.10f +/- %.10f E = %.10f ET = %.10f <E>_exact = %.10f \
                %6d Walkers \n", idmc, τ, Eb_ave, Eb_std, ETb_ave, ETb_std, E/n,  ET, Eexact, Nw)
            
            open(Efile,"a") do f
                println(f,τ," ",Eb_ave," ",Eb_std," ",Eexact," ",Nw)
            end
            
            if Nb>10 && Eb_std < accuracy_goal
                println("reached accuracy goal $accuracy_goal")
                println("result: τ  <E>  error")
                @printf("%g ",τ)
                open(Efile,"a") do f
                    println(f,τ," ",Eb_ave," ",Eb_std," ",Eexact)
                end
                output_MCresult(Eb_ave, Eb_std)
                break
            end
            GC.gc() # manual carbage collection
        end
        
        # 
        if idmc==0            
            println("THERMALIZATION ENDS")
            #Estat = init_stat(1, blocksize) # only needed if measuring during thermalization
            E_ave = E_ave/nE
            nE = 1
        end
    end
    
end

main()
