using Printf
#using Distributions
using StaticArrays
using Statistics
using ForwardDiff: gradient as ADgrad 
#using ReverseDiff: gradient as ADgrad 
#import LinearAlgebra: norm
# dot product:
import LinearAlgebra: ⋅  

# slightly faster (?) norm 
@inline norm(x) = sqrt(sum(abs2,x))

# local modules:
push!(LOAD_PATH,".")
using Common
using Jastrow: φ_T_J as φ_T_J_in, F_J as F_J_in
using VMCstep
using Utilities
using QMC_Statistics
using SlaterDeterminants
using LinearOptimization
using STOOrbitals

# for testing purposes only:
#using Random
#Random.seed!(1234)



possible_args=["atom","basis_set"]
arg_dict = argparse(possible_args)
# use command line values or defaults
atom = get(arg_dict, "atom", "")
basis_set = get(arg_dict, "basis_set", "STO")


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

jpara = parse_Jastrow_data(atom, basis_set)

println("Jastrow parameters:")
@show jpara.param

#sanity check
if jpara.param[4] == 0.0
    println("Warning: You have α12=0.0")
end
if jpara.param[2] == 0.0 || jpara.param[3] == 0.0
    println("Warning: You have β1 or β2 equal to 0.0")
end



# Select atomic orbital functions
SlaterDeterminants.init(basis_set)

all_aos = parse_basis_set_data(atom, basis_set)



# Read lc coefficients and nl-combinations from file  
cnls = parse_coeffients_data(atom, basis_set)

# create linear combination of ao's 
lc, npara_, nonlin_, lc_ups_, lc_dos_ = set_linear_combination_of_orbitals(atom_data, all_aos, cnls)

# fix lists of up and down orbitals
const lc_ups = lc_ups_
const lc_dos = lc_dos_
# sanity check
for aos in lc.orbitals
    if length(aos) != N
        @show length(aos),N
        error("lc.orbitals must have length N to accommondate all electrons")
    end
end


# Jastrow parameters
npara_ += 4  # hard-wired ...
const npara = npara_

const nonlin = nonlin_   # If you want to set this empty, use Int64[] 

optfile = string("atom_data/"*atom*"_"*basis_set*"_opt")
println("Optimized values will be saved to $optfile")


suffix_ = string(PROGRAM_FILE,"_",atom,"_opt")

const suffix = suffix_ ::String
const Efile = string("E.",suffix)
const VMCfile = string("VMC.",suffix)


println("="^40)
@show(suffix)
println("Energy output file: ",Efile)
println("   VMC output file: ",VMCfile)
println("="^40)

rm(Efile, force=true)
rm(VMCfile, force=true)

# bind Jastrow functions
φ_T_J(x ::MMatrix{dim, N, Float64}, par) = φ_T_J_in(x, par, s)
F_J(x ::MMatrix{dim, N, Float64}, par) = F_J_in(x, par, s)


# 2 for parallel spins, 4 for antiparallel
const s = spinfactor(atom_data.Nup, atom_data.Ndo)


# collect parameters to one struct
wf_params = WfParams(lc, jpara)
# erase temporaries
lc = nothing
jpara = nothing


println(npara," parameters in φ_T")
println(length(nonlin)," nonlinear parameters in φ_T")



@inline function E_correlated(walker ::Vector{Walker}, φ_T ::Function, wf_params_new::WfParams)
    E = 0.0        
    r = 0.0
    for w ∈ walker
        wratio = φ_T(w.R, wf_params_new)^2 / w.Ψ^2
        E += EL(w.R, wf_params_new) * wratio
        r += wratio
    end
    E = E/r
end



@inline function EL(R ::MMatrix{dim, N, Float64}, wf_params::WfParams)
    return EL(R, wf_params.lc, wf_params.jpara)
end

@inline function EL(R ::MMatrix{dim, N, Float64}, lc ::LinearCombination, jpara ::JastrowParameters)
    
    VL = 0.0 

    ∇DD = Matrix{Float64}(undef,dim,N)
    for j in 1:N
        for k in 1:dim
            ∇DD[k,j] = 0.0
        end
    end
    TL = 0.0
    DDsum = 0.0
    for (i,aos) in enumerate(lc.orbitals)        
        ups = aos[lc_ups[i]]
        dos = aos[lc_dos[i]]       
        # see comments before F(R)           
        DD_occu = lc.coeffs[i] * D_up(R, ups) * D_do(R, dos)
        TL += -D* DD_occu*( ∇2_up(R, ups) + ∇2_do(R, dos) )
        ∇DD +=   DD_occu *( F_up(R, ups) + F_do(R, dos) )
        DDsum +=  DD_occu
    end
    TL = TL/DDsum
    ∇DD = ∇DD/(DDsum*2)
    
    ∇J   = F_J(R, jpara)/2     # ∇ e^J/J = ∇J = F/2
    
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

#
# The e-n cusp condition is 
# c1*ζ1 + c2*ζ2 + ... = (Z-α)(c1+c2+...)
# with Jastrow factor parameter α
# For a single 1s orbital this reduces to ζ1 = Z.
# 

function force_en_cusp_condition!(wf_params ::WfParams)
    cs = []
    ζs = []
    for aos in wf_params.lc.orbitals
        done = false
        for ao in aos
            if [ao.n, ao.l] == [1, 0] 
                # this is 1s orbital
                cs = ao.coefficients
                ζs = ao.exponents
                done = true
                break
            end                
        end
        if done
            break
        end
    end
    α = wf_params.jpara.param[1]
    # sum(cs.*ζs) = (Z-α)*sum(cs), but not very accurate for many-electron systems

    #cs[end] = ((Z-α)*sum(cs[1:end-1]) - sum(cs[1:end-1].*ζs[1:end-1]) )/(ζs[end]-(Z-α))
    
    x  =  (Z-α)*sum(cs)/sum(cs.*ζs)
    println("cusp condition factor $x")
    if x < 0.0
        println("Warning: failed to enforce e-n cusp condition, would lead to ζ<0")
        @show  sum(cs.*ζs) - (Z-α)*sum(cs)
        return
    end
    ζs .*= x
    orbitals = []
    for aos in wf_params.lc.orbitals
        new_aolist = AtomicOrbital[]                
        for ao in aos
            if [ao.n, ao.l] == [1, 0]
                ao.exponents .= ζs
                ao.coefficients .= cs
            end
            push!(new_aolist, ao)
        end
        push!(orbitals, new_aolist)
    end
    wf_params.lc.orbitals .= orbitals
    
end

@inline function φ_T(R ::MMatrix{dim, N, Float64}, wf_params ::WfParams)
    return φ_T(R, wf_params.lc, wf_params.jpara)
end

# total Slater-Jastrow trial wavefunction

@inline function φ_T(R ::MMatrix{dim, N, Float64}, lc ::LinearCombination, jpara ::JastrowParameters) 
    DD = 0.0    
    for (i,aos) in enumerate(lc.orbitals)
        ups = aos[lc_ups[i]]
        dos = aos[lc_dos[i]]
        DD += lc.coeffs[i] * D_up(R, ups) * D_do(R, dos)
    end
    if isapprox(DD,0.0)
        error("φ_T: DD = 0, probably bad orbitals")
    end
    φ_T = φ_T_J(R, jpara) * DD
end


# version with default parameters
@inline @fastmath function φ_T(R ::MMatrix{dim, N, Float64}) 
    return φ_T(R, wf_params)
end

function G_i_numerical(R ::MMatrix{dim, N, Float64}, wf_params ::WfParams, G :: Function)
    G0 = G(R, wf_params)
    h = 1e-6
    G_i = zeros(npara)
    Δα = @MVector zeros(npara)
    for i in 1:npara
        Δα[i] = h
        wf_params = update_wf_params(wf_params, Δα)
        G_i[i] =  (G(R, wf_params) - G0)/h
        wf_params = update_wf_params(wf_params, -Δα)
        Δα .= 0.0
    end
    G_i  
end


# AD for both Ψ_i and EL_i computation; just different input G  
@inline function G_i_AD(R ::MMatrix{dim, N, Float64}, wf_params ::WfParams, G ::Function)

    G_J(x) = G(R, wf_params.lc, x)
    G_S(x) = G(R, x, wf_params.jpara)
   
    
    # Gradient wrt. lc.coeffs and parameters in each AtomicOrbital 
    function gradient_lc(lc)
        para_concat = Float64[]
        append!(para_concat, lc.coeffs)  # linear expansion coefficients
        
        nl = []
        for aos in lc.orbitals
            for ao in aos
                if [ao.n, ao.l] ∉ nl
                    push!(nl,[ao.n, ao.l])
                    append!(para_concat, ao.exponents)
                    append!(para_concat, ao.coefficients)
                end
            end
        end
        
        gradients = ADgrad(para_concat) do x
            pos = 1
            len = length(lc.coeffs)
            coeffs = x[pos:pos+len-1]
            pos += len
            
            nl = []
            new_aolist = AtomicOrbital[]                
            for aos in lc.orbitals                
                for ao in aos
                    if [ao.n, ao.l] ∉ nl
                        push!(nl,[ao.n, ao.l])
                        len = length(ao.exponents)    # coefficient have same len                    
                        new_exponents = x[pos:pos+len-1]
                        pos += len
                        new_coefficients = x[pos:pos+len-1]
                        pos += len
                        new_ao = AtomicOrbital(ao.atom, ao.orbital_type, ao.n, ao.l, new_exponents, new_coefficients, ao.spin)
                        push!(new_aolist, new_ao)
                    end
                end
            end
            # lc may have same ao multiple times, and there are fewer parameters to derive
            # Go through all orbitals in lc and pick their derived versions from  new_aolist
            orbitals = []
            for aos in lc.orbitals
                new_aos = AtomicOrbital[]      
                for ao in aos
                    # find the derived ao
                    for (this_nl, this_ao) in zip(nl,new_aolist)                        
                        if ao.n == this_nl[1] && ao.l == this_nl[2]
                            # pick *only* parameters from this_ao! (spin may differ)                        
                            new_ao = AtomicOrbital(ao.atom, ao.orbital_type, ao.n, ao.l, this_ao.exponents, this_ao.coefficients, ao.spin)
                            push!(new_aos, new_ao)
                            break
                        end
                    end                        
                end
                push!(orbitals, new_aos)
            end
            new_lc = LinearCombination(orbitals, coeffs)
            return G_S(new_lc)
        end
        return gradients
    end

    # wrapper function to differentiate wrt. struct fields
    function G_params(param)
        new_jpara = JastrowParameters(param)
        return G_J(new_jpara)
    end
    # Gradient wrt. jpara
    function gradient_jpara(jpara::JastrowParameters)
        param = jpara.param
        gradients = ADgrad(G_params, param)
        return gradients
    end
    
    grad_lc = gradient_lc(wf_params.lc)   
    grad_jpara  = gradient_jpara(wf_params.jpara)
    
    G_i = reduce(vcat, [grad_lc, grad_jpara]) # keep order!
    G_i
end


function save_parameters(outfile ::String, wf_params ::WfParams)
    #println("Saving parameters to file ",outfile)
    open(outfile,"w") do f
        println(f,atom*" "*basis_set)
        println(f,"*")
        order = ["1s","2s","2p","3s","3p","4s","3d","4p","5s","4d","5p","6s","4f","5d","6p"]

        
        for orbital in order
            found = false
            for aos in wf_params.lc.orbitals
                for ao in aos
                    if occursin(orbital, ao.orbital_type) # covers 2px etc.
                        mytype = ao.orbital_type[2]
                        len = length(ao.exponents)
                        println(f, len,"  ",mytype)
                        for i in 1:len
                            @printf(f,"%15.10f %15.10f\n", ao.exponents[i], ao.coefficients[i])
                        end
                        found = true
                        break
                    end
                    if found
                        break
                    end
                end
                if found
                    break
                end
            end            
        end
        println(f,"*")
        println(f,"# Jastrow params")        
        @printf(f,"%s %15.10f\n","α = ",wf_params.jpara.param[1])
        @printf(f,"%s %15.10f\n","β1 = ",wf_params.jpara.param[2])
        @printf(f,"%s %15.10f\n","β2 = ",wf_params.jpara.param[3])
        @printf(f,"%s %15.10f\n","α12 = ",wf_params.jpara.param[4])

        println(f,"# lc coefficients")        
        for (cc, aos) in zip(wf_params.lc.coeffs, wf_params.lc.orbitals)
            @printf(f,"coeff = %15.10f \n ",cc)
            @printf(f,"orbitals = ")
            for ao in aos
                @printf(f,"%d%d ",ao.n, ao.l)
            end
            @printf(f,"\n")
        end

        
        # raw print all wf_params
        println(f,"# raw output of wf_params")
        for cc in wf_params.lc.coeffs
            @printf(f,"%8.5f ", cc)
        end
        nl = []
        for aos in wf_params.lc.orbitals
            for ao in aos
                if [ao.n, ao.l] ∉ nl
                    push!(nl,[ao.n, ao.l]) 
                    for ζ in ao.exponents
                        @printf(f,"%8.5f ",ζ)
                    end
                    for c in ao.coefficients
                        @printf(f,"%8.5f ",c)
                    end
                end
            end
        end
        for p in wf_params.jpara.param
            @printf(f,"%8.5f ",p)
        end
        println(f,"")
    end
end

function thermalize(walker ::Vector{Walker}, nsteps ::Int64, vmc_params ::VMC_Params)
    # thermalize
    EL_ave = 0.0
    EL2_ave = 0.0
  
    Nw = length(walker)
    @inbounds for iw in 1:Nw 
        w = walker[iw]
        @inbounds for i in 1:nsteps
            # make sure walkers have Ψ set, or LinearOptimization fails
            w.Ψ = vmc_step!(w.R, vmc_params, (x->φ_T(x, wf_params)))
            if i%10 == 0
                adjust_step!(vmc_params)
            end
        end
        # make sure walkers have E set
        E = EL(w.R, wf_params)
        w.E = E
        EL_ave += E            
        EL2_ave += E^2   
    end        
    EL_ave /= Nw
    EL2_ave /= Nw
    EL_sigma = sqrt(abs(EL2_ave - EL_ave^2))

    return walker, EL_ave, EL_sigma 
end

function solve_Δα_wf_params!(walker ::Vector{Walker},
                             Ψ_is ::Matrix{Float64},
                             EL_is ::Matrix{Float64},
                             ELs ::Vector{Float64},                             
                             nonlin ::Vector{Int64},
                             fix_fact ::MVector,
                             a_opt ::Float64,
                             Δα ::MVector,
                             wf_params ::WfParams)

    Nw = length(walker)
    npara = length(Δα)
    wf_params_orig = copy(wf_params)
    wf_params_tmp  = copy(wf_params)
    while true
        # ************************************************************
        solve_Δα!(walker, Ψ_is, EL_is, ELs, a_opt, npara, nonlin, Δα)
        # ************************************************************
        # Add a small random variation
        #Δα .+= (rand(npara).-0.5)./(Nw^2) 
        # possibly fix some parameters
        Δα .*= fix_fact

       
        
        mα = maximum(abs.(Δα))
        if mα>1e-1
            # too large step
            a_opt *= 2.0
            continue
        end        
        #print_params(Δα,"Δα")        
        wf_params_tmp = update_wf_params(wf_params_orig, Δα)
        
        #force_en_cusp_condition!(wf_params_tmp)

        wf_params_tmp = norm_params!(wf_params_tmp, basis_set, Z)

        #print_params(wf_params_tmp,"wf_params_tmp")

        ok = true
        @inbounds for iw in 1:Nw
            if φ_T(walker[iw].R, wf_params_tmp) < 0.0
                # BAD
                println("φ_T changed sign (nodal cell changed or negative boson wf)")
                # undo update
                wf_params_tmp = copy(wf_params_orig)
                a_opt *=  2.0                            
                ok = false
                break
            end                        
        end
        if !ok; continue; end
        # all well, can use this 
        break                
    end
    return a_opt, Δα, wf_params_tmp
end

function run_VMC!(walker ::Vector{Walker}, vmc_params ::VMC_Params) ::Tuple{Float64, Float64}
    Nw = length(walker)
    println("thermalizing $Nw walkers using VMC")
    @inbounds for iw in 1:Nw
        R = walker[iw].R
        @inbounds for i in 1:100
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
    while true
        E_ave = 0.0 
        @inbounds for iw in 1:Nw 
            R = walker[iw].R
            vmc_step!(R, vmc_params, φ_T)
            E = EL(R, wf_params) 
            E_ave += E            
        end

        add_sample!(Estat, E_ave/Nw)
        adjust_step!(vmc_params)
        if Estat.finished            
            Eb_ave, Eb_std, E_inputvar2, Nb = get_stats(Estat)            
            @printf("VMC <E> = %.10f +/- %.10f \n", Eb_ave, Eb_std)
            
            
            
            # exit VMC
            if (Eb_std<0.001 && Nb>5) || Nb==20
                output_MCresult(Eb_ave, Eb_std)
                open(VMCfile,"a") do f
                    @printf(f,"%.10f  %.10f \n", Eb_ave, Eb_std)
                end
                open(optfile,"a") do f
                    output_MCresult(f, Eb_ave, Eb_std)
                end
                return Eb_ave, Eb_std
            end
            
            
        end        
    end
end

function main()
    #
    # Main program 
    # 
    # 
    vmc_params = VMC_Params(0,0,2.0)

    
    init_this(Nw,Nwx) = init_walkers(Nw, Nwx, dim, N; φ_T=φ_T, scale=2.5) 
    vmc_step_this!(R) = vmc_step!(R, vmc_params, φ_T)


    # -----------------
    # Check EL 
    Nw = 5
    print_params(wf_params,"wf_params")
    walker = init_this(Nw,Nw)
    println("Checking numerically EL against φ_T")
    
    for i in 1:Nw
        R = walker[i].R
        E = EL(R, wf_params)
        num_check_EL(R, E, φ_T, x->V(x, N, Z))
    end
    println("Full EL check passed")    
    println("checks passed")
    # ---------------------
    

    #
    # Optimization
    #
    opt_logfile = "Eopt_"*atom*"_"*basis_set
    
    println("Begin Optimization, see file $opt_logfile")
    rm(opt_logfile, force=true)
    rm("Eopt_EL_"*atom*"_"*basis_set, force=true)

    
    Δα = @MVector zeros(npara)

    # Fix some parameters
    # set 0 factor for fixed parameters, 1 else
    fix_fact = @MVector ones(npara)

    
    # fix_fact[length(wf_params.lc.coeffs)+1] = 0.0  # fix 1s exponent
    #fix_fact[end-3] = 0.0 # fix Jastrow α
    fix_fact[end] = 0.0 # fix Jastrow α12
    
    if occursin("EMA",basis_set)
        # fix exponents for EMA
        fix_fact[nonlin] .= 0.0
    end

    Nw = 1000
    
    while true
        Nw = Int(floor(Nw*1.5))
        if Nw>500000
            break
        end
        # using @MMatrix for large Nw may cause ERROR: LoadError: syntax: expression too large
        Ψ_is  = zeros(npara, Nw)
        EL_is = zeros(npara, Nw)
        ELs   = zeros(Nw)
        
        println("="^80)
        @show(Nw)
        println("="^80)
        
        walker = init_walkers(Nw, Nw, dim, N; φ_T=(x->φ_T(x, wf_params)), scale=1.5)

        
        println("-"^80)
        
        # thermalize over given number of VMC steps
        walker, EL_ave_old, EL_sigma_old = thermalize(walker, 100, vmc_params)
        
        @printf("k = %4d  Nw = %d <EL> = %.10f σ(EL) = %.10f \n", -1, Nw, EL_ave_old, EL_sigma_old)           
        a_opt = 10.0
        print_params(wf_params,"wf_params")
        
        ntries = 0
        ndownhill = 0
        for k ∈ 1:20
            global wf_params
            
            for iw ∈ 1:Nw
                R = walker[iw].R
                Ψ_is[:,iw]  =  G_i_AD(R, wf_params, φ_T)              
                EL_is[:,iw] = G_i_AD(R, wf_params, EL)
                #Ψ_is[:,iw] = G_i_numerical(R, wf_params, φ_T)              
                #EL_is[:,iw]= G_i_numerical(R, wf_params, EL)
                ELs[iw] =  EL(R, wf_params)
                walker[iw].E = ELs[iw]
                walker[iw].Ψ = φ_T(R, wf_params)
            end

            # solve wf_params and correlated E for three choices of a_opt
            # -----------------------------------------------------------
            a_opts = [0.1*a_opt, a_opt, 100*a_opt]
            E_opts = zeros(3)
                       
                        
            for i_opt ∈ 1:3
                a_opt_tmp = a_opts[i_opt]
                a_opt_tmp, Δα, wf_params_tmp = solve_Δα_wf_params!(walker, Ψ_is, EL_is, ELs, nonlin, fix_fact, a_opt_tmp, Δα, wf_params)
                a_opts[i_opt] = a_opt_tmp # may have changed in solve
                E_opts[i_opt] = E_correlated(walker, φ_T, wf_params_tmp)
            end
            println("Correlated energies for three values of a_opt:")
            @show E_opts
            @show a_opts

            a_opt_ok = true
            
            # find a_opt that gives the lowest correlated E using parabolic fit 
            A = [a_opts[1]^2 a_opts[1] 1;
                 a_opts[2]^2 a_opts[2] 1;
                 a_opts[3]^2 a_opts[3] 1]                           
            # fit to ax^2+bx+c 
            a, b, c = A\E_opts # fails if there are <3 unique a_opts values!
            a_opt = -b/(2a)
            if a < 1e-14 || a_opt<1e-12  # maximum in parabola or just very small a_opt
                # use a_opts that gave lowest E_opts 
                imin = argmin(E_opts)
                a_opt = a_opts[imin]
            end
            @show a_opt
            
            
            # solve wf_params for a_opt found above
            # -------------------------------------
            wf_params_orig = copy(wf_params)

            a_opt, Δα, wf_params = solve_Δα_wf_params!(walker, Ψ_is, EL_is, ELs, nonlin, fix_fact, a_opt, Δα, wf_params)

            #
            
            print_params(Δα,"Δα")
            print_params(wf_params,"wf_params")
            # -------------------------------------
            
            mα = maximum(abs.(Δα))
            
            EL_ave = 0.0
            EL2_ave = 0.0           
            for w in walker
                E = EL(w.R, wf_params)
                EL_ave += E
                EL2_ave += E^2                
            end
            EL_ave /= Nw
            EL2_ave /= Nw
            EL_sigma = sqrt(abs(EL2_ave - EL_ave^2))
            #
            # 
            if EL_ave < atom_data.Eexact - 5.0
                println("getting too low EL, breaking")                
                break
            end
            ok = true
            
            ntries += 1
            EL_corr_ave = E_correlated(walker, φ_T, wf_params)
            #if EL_ave <= EL_ave_old + EL_sigma_old/sqrt(Nw)  # && EL_sigma < 2*EL_sigma_old
            # decision to accept/reject based on correlated energy
            if EL_corr_ave <= EL_ave_old  #+ EL_sigma_old/Nw
                # BETTER
                ndownhill += 1                
                EL_ave_old = EL_ave
                EL_sigma_old = EL_sigma
                
            else                
                # WORSE
                ok = false                
                wf_params = wf_params_orig
            end
            open(opt_logfile,"a") do f
                println(f,EL_ave, " ", EL_sigma, " ", EL_corr_ave)
            end
            
            
            @printf("***** Accepted  %.5f %% \n", ndownhill*100.0/ntries)

            @printf("k = %4d  Nw = %d <EL> = %.10f σ(EL) = %.10f max Δα = %.3e  a_opt = %.6f, ok = %3d\n",
                    k, Nw, EL_ave, EL_sigma, mα, a_opt, ok)

           
            
            if mα<1e-10
                break
            end
            # short thermalization; new configs and new EL_ave_old
            # may be good for small Nw
            #println("running short VMC using current wfparams")
            #walker, EL_ave_old, EL_sigma_old = thermalize(walker, 10, vmc_params)
            
        end
        
        GC.gc() # manual carbage collection
        save_parameters(optfile, wf_params)
        # this may be omitted, takes a while:
        if ndownhill>0
            Eb_av, Eb_std = run_VMC!(walker, vmc_params)
            if Eb_std<1e-10
                println("perfect trial wave function found")
                break
            end
        end
        
    end
    
    println("Optimization done")
    println("Results in file",optfile*"_tmp")
    println("Move it to ",optfile," if you're satisfied")

    
end

@time main()
