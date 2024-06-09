using LinearAlgebra: norm, I, cholesky, Symmetric
using Serialization: deserialize
import Statistics

include("tril_indices.jl")

function load_data()
    Zatoms = deserialize("Zatoms.dat")
    R_all = deserialize("R.dat")
    E_all = deserialize("E.dat")
    F_all = deserialize("F.dat")
    idxs_train = deserialize("idxs_train.dat")
    #
    Ndata = length(E_all)
    @assert Ndata == length(R_all)
    @assert Ndata == length(F_all)
    #
    Natoms = length(Zatoms)
    @assert Natoms == size(R_all[1],2)
    @assert Natoms == size(F_all[1],2)
    #
    return Zatoms, R_all, E_all, F_all, idxs_train
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


function assemble_Kmatrix!(K, jrow::Int64, Natoms, R_desc_v, R_d_desc_v)

    σ = 20     # kernel parameter, should the the same as task["sig"] ???
    sqrt5 = sqrt(5)
    mat52_base_div = 3*σ^4
    sig_pow2 = σ^2

    Rj_d_desc_full = uncompress_R_d(Natoms, R_d_desc_v[jrow])
    
    dim_i = 3*Natoms
    Kij = zeros(Float64, dim_i, dim_i)

    idx_col_start = (jrow-1)*dim_i + 1
    idx_col_stop = jrow*dim_i
    idx_cols = idx_col_start:idx_col_stop

    for irow in jrow:Ntrain

        idx_row_start = (irow-1)*dim_i + 1
        idx_row_stop = irow*dim_i
        idx_rows = idx_row_start:idx_row_stop

        diff_ab = R_desc_v[irow] - R_desc_v[jrow]
        norm_ab = sqrt5 * norm(diff_ab)
        mat52_base = exp(-norm_ab/σ) / mat52_base_div * 5
    
        res1 = diff_ab * mat52_base * 5
        res2 = Rj_d_desc_full * diff_ab # matmul
        diff_ab_outer = res1 * res2'  # matmul
    
        res4 = (sig_pow2 + σ * norm_ab) * mat52_base  # scalar
        diff_ab_outer .-= res4 * Rj_d_desc_full'   # matmul
        Ri_desc_full = uncompress_R_d(Natoms, R_d_desc_v[irow])

        @views Kij[:,:] .= Ri_desc_full * diff_ab_outer
        @views K[idx_rows,idx_cols] .= Kij[:,:]
        @views K[idx_cols,idx_rows] .= Kij[:,:]'
    end

    return

end



Zatoms, R_all, E_all, F_all, idxs_train = load_data()
Natoms = length(Zatoms)

Ntrain = length(idxs_train)
# Calculate all descriptors
R_desc_v = Vector{Vector{Float64}}(undef,Ntrain)
R_d_desc_v = Vector{Matrix{Float64}}(undef,Ntrain)
ip = 1
for i in idxs_train
    R_desc_v[ip], R_d_desc_v[ip] = calc_descriptor(Natoms, R_all[i])
    ip += 1
end

Nrows = Ntrain * 3 * Natoms
Ncols = Nrows
K = zeros(Float64, Nrows, Ncols)
for jrow in 1:Ntrain
    assemble_Kmatrix!(K, jrow, Natoms, R_desc_v, R_d_desc_v)
end

y = zeros(Float64, Ntrain*3*Natoms)
dim_i = 3*Natoms
ip = 1
for itrain in idxs_train
    idx_start = (ip-1)*dim_i + 1
    idx_stop = ip*dim_i
    @views y[idx_start:idx_stop] .= F_all[itrain][:]
    ip += 1
end
y_std = Statistics.std(y, corrected=false) # ddof=0 in np.std
y *= (1/y_std)

K .*= -1
λ = 1e-10  # regularization strength
K[:,:] .= K[:,:] + I*λ

K_chol_fact = cholesky(Symmetric(K))
α_F = K_chol_fact\y
α_F .*= -1

α_F_res = reshape(α_F, (3, Natoms, Ntrain))

desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64
idx_rows, idx_cols, idx_lin = tril_indices(Natoms)

R_d_desc_α = zeros(Float64, desc_dim, Ntrain)
dvji = zeros(Float64, 3)
for itrain in 1:Ntrain
    ip = 1
    for (i, j) in zip(idx_rows, idx_cols)
        @views dvji[:] = α_F_res[:,j,itrain] - α_F_res[:,i,itrain]
        R_d_desc_α[ip,itrain] = dot(R_d_desc_v[itrain][:,ip], dvji)
        ip += 1
    end
end

