function gen_basis(Nsites, Nup, Ndn)

    NupTotal = binomial(Nsites,Nup)
    NdnTotal = binomial(Nsites,Ndn)
    NupdnTotal = NupTotal * NdnTotal

    upStates = Int64[]
    for i in 0:(2^Nsites-1)
        bitStr = string(i, base=2) # convert to binary string
        n_one = count( ==('1'), bitStr ) # count occurence of one
        if n_one == Nup
            push!(upStates, i)
        end
    end

    dnStates = Int64[]
    for i in 0:(2^Nsites-1)
        bitStr = string(i, base=2)
        n_one = count( ==('1'), bitStr ) # count occurence of one
        if n_one == Ndn
            push!(dnStates, i)
        end
    end

    combinedBasis = zeros(Int64,NupdnTotal,2)
    i = 0
    # XXX This is the convention for the double-loop, may consider another choices
    for iup in 1:NupTotal, idn in 1:NdnTotal
        i = i + 1
        combinedBasis[i,1] = upStates[iup]
        combinedBasis[i,2] = dnStates[idn]
    end

    return upStates, dnStates, combinedBasis
end