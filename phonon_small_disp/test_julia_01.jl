using Infiltrator
using LinearAlgebra: diagm
import JSON

function read_forces_json(filename)
    f = open(filename, "r")
    line = readline(f)
    parsed_data = JSON.parse(line)["forces"]["__ndarray__"]
    array_shape = convert(Vector{Int64}, parsed_data[1])
    forces = convert(Vector{Float64}, parsed_data[3])
    forces = reshape(forces, (array_shape[1], array_shape[2]))
    return Matrix(forces') # transpose and convert to matrix
end

function gen_R_latvecs_idx(Ncell)
    NcellTotal = prod(Ncell)
    R_latvecs_idx = zeros(Int64, 3, NcellTotal)
    Nx = Ncell[1]
    Ny = Ncell[2]
    Nz = Ncell[3]
    ip = 1
    for i in 0:Nx-1, j in 0:Ny-1, k in 0:Nz-1
        R_latvecs_idx[:,ip] .= [i, j, k]
        ip += 1
    end
end

function main_debug()
    atsymbs = ["H", "H"]
    Natoms = 2
    atpos = zeros(Float64, 3, Natoms)
    atpos[:,1] = [0.0, 0.0, 0.0]
    atpos[:,2] = [0.91, 0.0, 0.0]
    unit_latvecs = diagm([1.8, 20.0, 20.0])
    Δ = 0.01
    idx_moved = [1, 2] # start from 1
    dirs = ["x", "y", "z"]
    Ncell = [9, 1, 1]
    file_prefix_name = "TEMP_phonon_H2_01"
    C_av_mat = Matrix{Matrix{Float64}}(undef, 3, Natoms)
    for (i,a) in enumerate(idx_moved)
        for (j,v) in enumerate(dirs)
            filename_m = joinpath(file_prefix_name, "cache."*string(a-1)*v*"-.json")
            filename_p = joinpath(file_prefix_name, "cache."*string(a-1)*v*"+.json")
            println("i=$i a=$a j=$(j) v=$(v) file=$(filename_m)")
            forces_m = read_forces_json(filename_m)
            forces_p = read_forces_json(filename_p)
            C_av_mat[j,i] = (forces_m - forces_p)/(2*Δ)
        end
    end
    @infiltrate
    return
end

#main_debug()

