include("dec2binstr_padded.jl")
include("gen_basis.jl")

function main()
    Nsites = 4
    Nup = 3
    Ndn = 2
    upStates, dnStates, combinedBases = gen_basis(Nsites, Nup, Ndn)
    println("upStates = ", upStates)
    println("dnStates = ", dnStates)
    println("combinedBases = ")
    display(combinedBases); println()
end

main()