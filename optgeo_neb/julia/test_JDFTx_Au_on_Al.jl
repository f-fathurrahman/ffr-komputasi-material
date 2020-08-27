using Printf
using LinearAlgebra

include("constants.jl")
include("Atoms.jl")
include("JDFTx.jl")

function main()

    atoms = Atoms(xyz_string="""
    13

    Al     0.00000000000000     0.00000000000000     4.00000000000000
    Al     2.86378246380552     0.00000000000000     4.00000000000000
    Al     0.00000000000000     2.86378246380552     4.00000000000000
    Al     2.86378246380552     2.86378246380552     4.00000000000000
    Al     1.43189123190276     1.43189123190276     6.02500000000000
    Al     4.29567369570828     1.43189123190276     6.02500000000000
    Al     1.43189123190276     4.29567369570828     6.02500000000000
    Al     4.29567369570828     4.29567369570828     6.02500000000000
    Al     0.00000000000000     0.00000000000000     8.05000000000000
    Al     2.86378246380552     0.00000000000000     8.05000000000000
    Al     0.00000000000000     2.86378246380552     8.05000000000000
    Al     2.86378246380552     2.86378246380552     8.05000000000000
    Au     1.43189123190276     1.43189123190276     9.75000000000000
    """,
    LatVecs=ANG2BOHR*diagm([5.72756492761104, 5.72756492761104, 13.75])
    )
    #println(atoms)

    calc = JDFTxCalculator()
    calc.Ncore = 2
    calc.use_smearing = true
    calc.kpoint_folding[:] = [3,3,1]
    calc.prefix_dir = "rundir_jdftx_Au_on_Al"
    forces = zeros(3,atoms.Natoms)
    energy = compute!(calc, atoms, forces)
    
    println("energy (Ha) = ", energy)
    println("forces = ")
    display(forces'); println()
end

main()