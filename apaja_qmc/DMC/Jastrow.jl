#
# Jastrow Factor 
#
module Jastrow

using StaticArrays

push!(LOAD_PATH,".")
using Common

export φ_T_J,  lnφ_T_J, F_J

# slightly faster norm 
@inline norm(x) = sqrt(sum(abs2,x))

# Jastrow trial wave function
@inline function φ_T_J(R ::MMatrix{dim,N,Float64}, jpara::JastrowParameters, s ::Matrix{Int64}) where {dim, N}
    return exp(lnφ_T_J(R, jpara, s))
end



@inline function lnφ_T_J(R ::MMatrix{dim,N,Float64},  jpara::JastrowParameters, s ::Matrix{Int64}) where {dim, N}
    lnφ_T = 0.0 
    α, β1, β2, α12 = jpara.param
    for k in 1:N
        rk  = norm(R[:,k])
        lnφ_T += -α*rk
        for i in k+1:N
            ski = s[k,i]
            rki = norm(R[:,k]-R[:,i])
            if ski==2
                lnφ_T += α12*rki/(ski*(1+β1*rki))
            else
                lnφ_T += α12*rki/(ski*(1+β2*rki))
            end
        end
    end
    return lnφ_T
end



# Jastrow drift
# F = 2∇J = 2∇_k φ_T /φ_T, Jastrow φ_T = e^J
@inline function F_J(R ::MMatrix{dim,N,Float64}, jpara::JastrowParameters, s ::Matrix{Int64})  where {dim, N}
    α, β1, β2, α12 = jpara.param
    ∇J = Array{eltype(jpara.param),2}(undef,dim,N)
    for i in 1:N
        for k in 1:dim
            ∇J[k,i] = 0.0
        end
    end
    
    for k ∈ 1:N
        rk = norm(R[:,k])
        ∇J[:,k] += -α*R[:,k]/rk
        for i ∈ 1:N
            if i==k; continue; end
            ski = s[k,i]
            rki_vec = (R[:,k]-R[:,i])
            rki = norm(rki_vec)
            if ski==2
                ∇J[:,k] +=  α12/(ski*(1+β1*rki)^2) * rki_vec/rki
            else
                ∇J[:,k] +=  α12/(ski*(1+β2*rki)^2) * rki_vec/rki
            end
        end
    end
    2∇J
end

end
