using LinearAlgebra: norm
using Serialization: deserialize

include("tril_indices.jl")


R_all = deserialize("R.dat")
E_all = deserialize("E.dat")
F_all = deserialize("F.dat")
Zatoms = deserialize("Zatoms.dat")

Ndata = length(E_all)
@assert Ndata == length(R_all)
@assert Ndata == length(F_all)

Natoms = length(Zatoms)
@assert Natoms == size(R_all[1],2)
@assert Natoms == size(F_all[1],2)

idata = 1
R = R_all[1]

dim_i = 3*Natoms
desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64

ip = 1
dR = zeros(Float64, desc_dim)
for icol in 1:Natoms, irow in (icol+1):Natoms
    dR[ip] = norm(R[:,irow] - R[:,icol])
    ip += 1 
end
R_desc = 1 ./ dR

R_d_desc = zeros(Float64, 3, desc_dim)
ip = 1
dR_vec = zeros(Float64, 3)
for icol in 1:Natoms, irow in (icol+1):Natoms
    @views dR_vec[:] .= R[:,irow] - R[:,icol]
    @views R_d_desc[:,ip] .= dR_vec[:] / dR[ip]^3
    ip += 1 
end
# original Python code have opposite sign

# Uncompressed R_d
function uncompress_R_d(R_d_desc)
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


