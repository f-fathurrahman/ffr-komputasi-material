using Printf

include("Atoms.jl")
include("MullerBrown.jl")

function main()
    mb = MullerBrown()
    forces = zeros(Float64,3,1)
    energy = calc_energy_forces!(mb, -0.55, 1.30, forces)
    println(energy)
    println(forces)
end

main()