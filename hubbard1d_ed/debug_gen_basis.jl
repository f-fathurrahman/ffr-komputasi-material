include("dec2binstr_padded.jl")
include("gen_basis.jl")

function main()
    Nsites = 4
    Nup = 3
    Ndn = 2
    upStates, dnStates, combinedBasis = gen_basis(Nsites, Nup, Ndn)
    println("upStates = ", upStates)
    println("dnStates = ", dnStates)
    println("combinedBasis = ")
    display(combinedBasis); println()
end

main()