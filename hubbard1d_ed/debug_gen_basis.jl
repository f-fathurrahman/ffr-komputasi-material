
function main()
    Nsites = 4
    Nup = 3
    Ndn = 2

    NupTotal = binomial(Nsites,Nup)
    NdnTotal = binomial(Nsites,Ndn)
    NupdnTotal = NupTotal * NdnTotal

    upStates = Int64[]
    for i in 0:(2^Nsites-1)
        bitStr = string(i, base=2) # convert to binary string
        n_one = count( ==('1'), bitStr ) # count occurence of one
        if n_one == Nup
            #println("bitStr = ", bitStr)
            #println("n_one = ", n_one, " i = ", i)
            push!(upStates, i)
        end
    end
    println("upStates = ", upStates)

    dnStates = Int64[]
    for i in 0:(2^Nsites-1)
        bitStr = string(i, base=2)
        n_one = count( ==('1'), bitStr ) # count occurence of one
        if n_one == Ndn
            #println("bitStr = ", bitStr)
            #println("n_one = ", n_one, " i = ", i)
            push!(dnStates, i)
        end
    end
    println("dnStates = ", dnStates)

    combinedBases = zeros(Int64,NupdnTotal,2)
    i = 0
    for iup in 1:NupTotal, idn in 1:NdnTotal  # This is
        i = i + 1
        combinedBases[i,1] = upStates[iup]
        combinedBases[i,2] = dnStates[idn]
    end
    display(combinedBases); println()

end

main()