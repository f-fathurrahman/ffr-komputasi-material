include("MullerBrown.jl")

function main()
    mb = MullerBrown()
    energy, forces = calc_energy_forces(mb, -0.55, 1.30)
    println(energy)
    println(forces)
end

main()