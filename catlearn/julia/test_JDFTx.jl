using Printf
using LinearAlgebra

include("constants.jl")
include("Atoms.jl")
include("JDFTx.jl")

function main()
    atoms = Atoms(xyz_string="""
    4

    O   0.0  1.0   0.0
    H   0.0  0.0   0.0
    H   1.1  0.0   0.0
    C   0.0  0.0   1.0
    """, LatVecs=16.0*diagm([1.0, 1.0, 1.0])
    )
    #println(atoms)

    calc = JDFTxCalculator()
    calc.Ncore = 2
    calc.use_smearing = true
    forces = zeros(3,atoms.Natoms)
    energy = compute!(calc, atoms, forces)
    
    println("energy (Ha) = ", energy)
    println("forces = ")
    display(forces'); println()
end

main()