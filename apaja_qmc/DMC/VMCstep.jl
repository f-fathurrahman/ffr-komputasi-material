module VMCstep

using StaticArrays
using Random: rand!

push!(LOAD_PATH,".")
using Common: VMC_Params
using Utilities: dist1, dist2
export vmc_step!, adjust_step!


const buf = Vector{Float64}(undef, 1024 * 3) # max_N * max_dim
const d_buf = MVector{3, Float64}(undef)


# alternatives
# ------------
function vmc_step!(R::MMatrix{dim, N, Float64}, params ::VMC_Params, Ψ::Function) where {dim,N}
    # sanity checks:
    @assert N <= 1024
    if params.step<1e-15
        error("VMC step is zero")
    end
    ΨR = Ψ(R)
    if ΨR < 0.0
        error("Ψ(R)<0, this shouldn't happen.")
    end
    Ψ2 = ΨR^2
    rrs = reshape(view(buf, 1:(dim*N)), dim, N)
    rand!(rrs)
    d = @view d_buf[1:dim]
    @inbounds for i in 1:N
        rr = @view rrs[:, i]
        @. d  = params.step*( rr - 0.5 )
        R[:, i] += d
        Ψ_new = Ψ(R)
        Ψ2_new = Ψ_new^2
        ratio = Ψ2_new/Ψ2
        # Metropolis:
        accept = false        
        params.Ntry +=1
        if Ψ_new<0.0
            # fermi nodal cell change, and shoudn't happen for boson ground state either
            accept = false
        else        
            if ratio >= 1.0
                accept = true            
            else
                if ratio > rand() 
                    accept = true
                end
            end
        end
        if accept
            params.Naccept +=1
            Ψ2 = Ψ2_new            
        else
            # revert to old value
            R[:, i] -= d
        end        
    end
    sqrt(Ψ2)
end
#
# Electrons and protons moving, first half electrons, second half protons
# Used in H2 code
# Not to be used with a wave function that has nodes
#
function vmc_step!(R ::MMatrix{dim,N}, params ::Vector{VMC_Params}, Ψ::Function) where {dim,N}
    ΨR = Ψ(R)
    Ψ2 = ΨR^2
    Nhalf = Int(N/2)
    rrs = reshape(view(buf, 1:(dim*N)), dim, N)
    rand!(rrs)
    d = @view d_buf[1:dim]
    # electrons and protons move with different steps
    for i in 1:N
        p = (i <= Nhalf) ? params[1] : params[2]
        rr = @view rrs[:, i]
        @. d  = p.step*( rr - 0.5 )
        R[:, i] += d
        Ψ_new = Ψ(R)
        Ψ2_new = Ψ_new^2
        ratio = Ψ2_new/Ψ2
        # Metropolis:
        accept = false        
        p.Ntry +=1        
        if ratio >= 1.0
            accept = true            
        elseif ratio > rand() 
            accept = true
        end
        if accept
            p.Naccept +=1
            Ψ2 = Ψ2_new            
        else
            # revert to old value
            R[:, i] -= d
        end        
    end
    sqrt(Ψ2)
end
#
# Only electrons moving
# Used in H2 code
# Not to be used with a wave function that has nodes
#

export vmc_step_H2!
    
function vmc_step_H2!(R ::MMatrix{dim, N}, p ::VMC_Params, Ψ::Function) where {dim, N}
    ΨR = Ψ(R)
    Ψ2 = ΨR^2
    Nhalf = Int(N/2)
    rrs = reshape(view(buf, 1:(dim*N)), dim, N)
    rand!(rrs)
    d = @view d_buf[1:dim]
    for i in 1:Nhalf
        rr = @view rrs[:, i]
        @. d  = p.step*( rr - 0.5 )
        R[:, i] += d
        Ψ_new = Ψ(R)
        Ψ2_new = Ψ_new^2
        ratio = Ψ2_new/Ψ2
        # Metropolis:
        accept = false        
        p.Ntry +=1        
        if ratio >= 1.0
            accept = true            
        elseif ratio > rand() 
            accept = true
        end
        if accept
            p.Naccept +=1
            Ψ2 = Ψ2_new            
        else
            # revert to old value
            R[:, i] -= d
        end        
    end
    sqrt(Ψ2)
end



function vmc_step!(R ::MMatrix{dim,N}, params ::VMC_Params, wf_params ::Vector{Float64}, Ψ::Function) where {dim,N}
    ΨR(x) = Ψ(x, wf_params)
    Ψ2 = ΨR(R)^2
    rrs = reshape(view(buf, 1:(dim*N)), dim, N)
    rand!(rrs)
    d = @view d_buf[1:dim]
    @inbounds for i in 1:N
        rr = @view rrs[:, i]
        @. d   = params.step*( rr - 0.5 )
        R[:, i] += d
        Ψ2_new = ΨR(R)^2
        ratio = Ψ2_new/Ψ2
        # Metropolis:
        accept = false        
        params.Ntry +=1
        if ratio >= 1.0
            accept = true            
        else
            if ratio > rand() 
                accept = true
            end
        end
        if accept
            params.Naccept +=1
            Ψ2 = Ψ2_new            
        else
            # revert to old value
            R[:, i] -= d
        end        
    end
    sqrt(Ψ2)
end




const minstep = 1e-5
const maxstep = 20.0

# adjust step to keep acceptance 50-60 %
function adjust_step!(params ::VMC_Params)
    acceptance = params.Naccept*100.0/params.Ntry
    if acceptance<50.0
        params.step *= 0.9
    end
    if acceptance>60.0
        params.step *= 1.1
    end
    params.step = max(minstep,params.step)
    params.step = min(maxstep,params.step)
    if params.step == maxstep
        println("BAD Error: step is maxstep ",maxstep)
        @show(acceptance)
        exit()
    end
end


end
