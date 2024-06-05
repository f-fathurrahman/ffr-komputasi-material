using LinearAlgebra: norm
using Serialization: deserialize

include("tril_indices.jl")

function load_data()
    Zatoms = deserialize("Zatoms.dat")
    R_all = deserialize("R.dat")
    E_all = deserialize("E.dat")
    F_all = deserialize("F.dat")
    #
    Ndata = length(E_all)
    @assert Ndata == length(R_all)
    @assert Ndata == length(F_all)
    #
    Natoms = length(Zatoms)
    @assert Natoms == size(R_all[1],2)
    @assert Natoms == size(F_all[1],2)
    #
    return Zatoms, R_all, E_all, F_all
end


# Uncompressed R_d
function uncompress_R_d(Natoms, R_d_desc)
    desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64
    tmp_R_d = zeros(Float64, 3, Natoms, desc_dim)
    idx_rows, idx_cols, idx_lin = tril_indices(Natoms)
    for ip in 1:desc_dim
        i = idx_rows[ip]
        j = idx_cols[ip]
        @views tmp_R_d[:,i,ip] .=  R_d_desc[:,ip]
        @views tmp_R_d[:,j,ip] .= -R_d_desc[:,ip]
    end
    R_d_full = reshape(tmp_R_d, 3*Natoms, desc_dim)
    # Use SparseArray ???
    return R_d_full
end


function calc_descriptor(Natoms, R)
    dim_i = 3*Natoms
    desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64
    #
    ip = 1
    dR = zeros(Float64, desc_dim)
    for icol in 1:Natoms, irow in (icol+1):Natoms
        dR[ip] = norm(R[:,irow] - R[:,icol])
        ip += 1 
    end
    R_desc = 1 ./ dR
    #
    R_d_desc = zeros(Float64, 3, desc_dim)
    ip = 1
    dR_vec = zeros(Float64, 3)
    for icol in 1:Natoms, irow in (icol+1):Natoms
        @views dR_vec[:] .= R[:,irow] - R[:,icol]
        @views R_d_desc[:,ip] .= dR_vec[:] / dR[ip]^3
        ip += 1 
    end
    # original Python code have opposite sign
    return R_desc, R_d_desc
end


Zatoms, R_all, E_all, F_all = load_data()
Natoms = length(Zatoms)

Ntrain = 2
R_desc_v = Vector{Vector{Float64}}(undef,Ntrain)
R_d_desc_v = Vector{Matrix{Float64}}(undef,Ntrain)

for i in 1:2
    R_desc_v[i], R_d_desc_v[i] = calc_descriptor(Natoms, R_all[i])
end


σ = 20     # kernel parameter, should the the same as task["sig"] ???
λ = 1e-10  # regularization strength, should be the same as task["lam"] ???

using LinearAlgebra: norm

i = 2
j = 2

sqrt5 = sqrt(5)
mat52_base_div = 3*σ^4
sig_pow2 = σ^2

diff_ab = R_desc_v[i] - R_desc_v[j]
norm_ab = sqrt5 * norm(diff_ab)
mat52_base = exp(-norm_ab/σ) / mat52_base_div * 5

res1 = diff_ab * mat52_base * 5
Rj_d_desc_full = uncompress_R_d(Natoms, R_d_desc_v[j])
res2 = Rj_d_desc_full * diff_ab # matmul
diff_ab_outer = res1 * res2'  # matmul

res4 = (sig_pow2 + σ * norm_ab) * mat52_base  # scalar
diff_ab_outer .-= res4 * Rj_d_desc_full'   # matmul

Ri_desc_full = uncompress_R_d(Natoms, R_d_desc_v[i])
Kij = Ri_desc_full * diff_ab_outer   # matmul
println("sum abs K = ", sum(abs.(Kij)))

