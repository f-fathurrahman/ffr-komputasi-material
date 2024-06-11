using Printf
using LinearAlgebra: norm, I, cholesky, Symmetric, dot
using Serialization: deserialize
import Statistics

include("tril_indices.jl")
include("GDMLModel.jl")

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

include("calc_descriptor.jl")
include("assemble_Kmatrix.jl")


function init_GDML_model()

    # Load data
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

    λ = 1e-10  # regularization strength
    σ = 20.0   # kernel parameter

    c = 0.0 # integration constant

    desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64
    R_d_desc_α = zeros(Float64, desc_dim, Ntrain)

    return GDMLModel(
        Natoms, Ntrain, desc_dim,
        R_desc_v,
        R_d_desc_v,
        σ, λ,
        R_d_desc_α
        TrilIndices(Natoms),
        c,
        y_std
    )


end



    Nrows = Ntrain * 3 * Natoms
    Ncols = Nrows
    K = zeros(Float64, Nrows, Ncols)
    for jrow in 1:Ntrain
        assemble_Kmatrix!(K, jrow, Natoms, R_desc_v, R_d_desc_v)
    end

    K .*= -1
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
        for (ia, ja) in zip(idx_rows, idx_cols)
            @views dvji[:] = α_F_res[:,ja,itrain] - α_F_res[:,ia,itrain]
            R_d_desc_α[ip,itrain] = dot(R_d_desc_v[itrain][:,ip], dvji)
            ip += 1
        end
    end

    return
end



function predict_train(itrain)

    # Predict for one data point
    r = R_all[idxs_train[itrain]] # just some some data from R_all
    #
    r_desc, r_d_desc = calc_descriptor(Natoms, r)

    diff_ab = zeros(Float64, desc_dim, Ntrain)
    norm_ab = zeros(Float64, Ntrain)
    for itrain in 1:Ntrain
        @views diff_ab[:,itrain] .= r_desc .- R_desc_v[itrain]
        norm_ab[itrain] = sqrt(5) * norm(diff_ab[:,itrain])
    end

    mat52_base = 5.0/(3*σ^3) * exp.(-norm_ab/σ)

    a_x2 = zeros(Float64, Ntrain)
    for itrain in 1:Ntrain
        @views a_x2[itrain] = dot(diff_ab[:,itrain], R_d_desc_α[:,itrain])
    end

    ff = diff_ab * (a_x2 .* mat52_base) * 5 / σ
    mat52_base .*= norm_ab .+ σ

    ff .-= R_d_desc_α * mat52_base

    E_pred0 = dot(a_x2, mat52_base)*y_std
    println("E_pred0 = ", E_pred0)

    # Here r_d_desc is used
    out_F = zeros(Float64, 3, Natoms, Natoms)
    ip = 1
    for (ia, ja) in zip(idx_rows, idx_cols)
        for i in 1:3
            out_F[i,ia,ja] = r_d_desc[i,ip] * ff[ip]
            out_F[i,ja,ia] = -out_F[i,ia,ja]
        end
        ip += 1
    end
    F_pred = dropdims(sum(out_F, dims=2), dims=2) * y_std
    # We sum over 2nd dimension here (to get the same sign for forces)
    println("F_pred = ")
    for ia in 1:Natoms
        @printf("%18.10f %18.10f %18.10f\n", F_pred[1,ia], F_pred[2,ia], F_pred[3,ia])
    end

    return
end

main()

