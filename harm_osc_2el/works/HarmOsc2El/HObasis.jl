struct HOBasis <: SpatialBasis
    l::Int64  # number of basis functions
    ω::Float64  # strength of harmonic oscillator potential
    hermites::Vector{Float64}
    hos::Vector{Float64}
    ho_der::Vector{Float64}
    ho_dder::Vector{Float64}

    function HOBasis(l, ω)
        hermites = zeros(l)
        hermites[1] = 1.0
        
        return new(l, ω, hermites, zeros(l), zeros(l), zeros(l))
    end
end

function evaluate(x, ho::HOBasis)
    (; ω, hermites) = ho
    hos = zero(hermites)
    
    x = √ω * x

    hermites[1] = 1.0
    hermites[2] = 2x

    ho_fac = (ω / π)^0.25 * exp(-x^2 / 2)
    hos[1] = ho_fac * hermites[1]
    ho_fac *= 1 / √2
    hos[2] = ho_fac * hermites[2]

    @inbounds for n in 3:length(hos)
        hermites[n] = 2x * hermites[n-1] - 2(n - 2) * hermites[n-2]

        ho_fac *= 1 / sqrt( 2(n - 1) )
        hos[n] = ho_fac * hermites[n]
    end
    
    return hos
end

function evaluate!(hos, x, ho::HOBasis)
    (; ω, hermites) = ho
    
    x = √ω * x

    #hermites[1] = 1.0
    ho_fac = (ω / π)^0.25 * exp(-x^2 / 2)

    hos[1] = ho_fac * hermites[1]

    hermites[2] = 2x
    ho_fac *= 1 / √2

    hos[2] = ho_fac * hermites[2]

    @inbounds for n in 3:length(hos)
        hermites[n] = 2x * hermites[n-1] - 2(n - 2) * hermites[n-2]
        ho_fac *= 1 / sqrt( 2(n - 1) )

        hos[n] = ho_fac * hermites[n]
    end
    
    return hos
end

function spatial(ho::HOBasis, grid)
    (; l, hos) = ho
    n = length(grid)
    res = [zeros(n) for i in 1:l]
    
    for i in 1:n
        evaluate!(hos, grid[i], ho)
        for j in 1:l
            @inbounds res[j][i] = hos[j]
        end
    end
    
    #for i in 1:l # In case I need the same spfs as ODQD from quantum_systems
    #    if i % 4 == 0 || (i + 1)%4 == 0
    #        res[i] = -res[i]
    #    end
    #end
    
    return res
end

function fast_ho_all!(x, ho::HOBasis)
    (; ω, hermites, hos, ho_der, ho_dder) = ho
    
    ξ = √ω * x

    #hermites[1] = 1.0
    ho_fac = (ω / π)^0.25 * exp(-ξ^2 / 2)
    
    hos[1] =     ho_fac * hermites[1]
    ho_der[1] = -ω * x * hos[1]
    ho_dder[1] = ω * (ω * x^2 - 1) * hos[1]
    
    hermites[2] = 2ξ
    ho_fac *= 1 / √2

    hos[2] =     ho_fac * hermites[2]
    ho_der[2] =  ho_fac * (√ω * 2 * hermites[1] - ω * x * hermites[2])
    ho_dder[2] = ho_fac * ω * ((ω * x^2 - 1) * hermites[2] - √ω * x * 4 * hermites[1])

    @inbounds for n in 2:length(hos)-1
        hermites[n+1] = 2ξ * hermites[n] - 2(n-1) * hermites[n-1]
        ho_fac *= 1 / sqrt( 2n )

        hos[n+1] = ho_fac * hermites[n+1]
        ho_der[n+1] =  ho_fac * (√ω * 2 * n * hermites[n] - ω * x * hermites[n+1])
        ho_dder[n+1] = ho_fac * ω * ((ω * x^2 - 1) * hermites[n+1] - √ω * x * 4n * hermites[n] + 4 * (n-1)*n * hermites[n - 1])
    end
    
    return hos, ho_der, ho_dder
end

function fast_ho_all!(hermites, hos, ho_der, ho_dder, x, ho::HOBasis)
    (; ω) = ho
    
    ξ = √ω * x

    hermites[1] = 1.0
    ho_fac = (ω / π)^0.25 * exp(-ξ^2 / 2)
    
    hos[1] =     ho_fac * hermites[1]
    ho_der[1] = -ω * x * hos[1]
    ho_dder[1] = ω * (ω * x^2 - 1) * hos[1]
    
    hermites[2] = 2ξ
    ho_fac *= 1 / √2

    hos[2] =     ho_fac * hermites[2]
    ho_der[2] =  ho_fac * (√ω * 2 * hermites[1] - ω * x * hermites[2])
    ho_dder[2] = ho_fac * ω * ((ω * x^2 - 1) * hermites[2] - √ω * x * 4 * hermites[1])

    @inbounds for n in 2:length(hos)-1
        hermites[n+1] = 2ξ * hermites[n] - 2(n-1) * hermites[n-1]
        ho_fac *= 1 / sqrt( 2n )

        hos[n+1] = ho_fac * hermites[n+1]
        ho_der[n+1] =  ho_fac * (√ω * 2 * n * hermites[n] - ω * x * hermites[n+1])
        ho_dder[n+1] = ho_fac * ω * ((ω * x^2 - 1) * hermites[n+1] - √ω * x * 4n * hermites[n] + 4 * (n-1)*n * hermites[n - 1])
    end
    
    return hos, ho_der, ho_dder
end