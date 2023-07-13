using Printf

include("gen_basis.jl")
include("dec2binstr_padded.jl")

function debug_kinetic()
    Nsites = 4
    Nup = 3
    Ndn = 2
    upStates, dnStates, combinedBasis = gen_basis(Nsites, Nup, Ndn)

    totalNoOfPossiblestates = size(combinedBasis,1)
    for m in 1:totalNoOfPossiblestates
        # save the unshifted spin up and spin down sectors:
        upSectorDec = combinedBasis[m,1]
        dnSectorDec = combinedBasis[m,2]

        #upSector = dec2binstr_padded(upSectorDec, Nsites)
        #dnSector = dec2binstr_padded(dnSectorDec, Nsites)
        #upNonZero = count(==('1'), upSector)
        #dnNonZero = count(==('1'), dnSector)

        upSector = reverse(digits(upSectorDec, base=2, pad=Nsites))
        dnSector = reverse(digits(dnSectorDec, base=2, pad=Nsites))

        upNonZero = findall(==(1), upSector)
        dnNonZero = findall(==(1), dnSector)

        # find the occupied lattice sites
        println("upNonZero = ", upNonZero)
        println("dnNonZero = ", dnNonZero)

        # shift for spin up:
        for n in upNonZero # for each occupied site
            # left shift:
            leftShiftResult = copy(upSector)
            # figure out which site is the one to its left (periodic boundary condition)
            leftShiftedIndex = mod(n-2,Nsites) + 1
            @printf("m = %3d n = %3d\n", m, n)
            #if upSector(leftShiftedIndex) ~= 1
        end

    end
end

debug_kinetic()
