using NPZ

filepath = "../gdml_python/DATASET/ethanol_dft.npz"

F_read = npzread(filepath, ["F"])["F"]
R_read = npzread(filepath, ["R"])["R"]
E_read = npzread(filepath, ["E"])["E"]
# Convert to Int64
Zatoms = Int64.( npzread(filepath, ["z"])["z"] )

Ndata = size(E_read, 1)
Natoms = size(R_read, 2)
@assert Natoms == size(F_read, 2)

E = dropdims(E_read, dims=2) # remove 2nd dimension (it is a singleton dimension)
R = Vector{Matrix{Float64}}(undef,Ndata)
F = Vector{Matrix{Float64}}(undef,Ndata)
for idata in 1:Ndata
    R[idata] = zeros(Float64,3,Natoms)
    F[idata] = zeros(Float64,3,Natoms)
    for ia in 1:Natoms, i in 1:3
        R[idata][i,ia] = R_read[idata,ia,i]
        F[idata][i,ia] = F_read[idata,ia,i]
    end
end

using Serialization
serialize("E.dat", E)
serialize("R.dat", R)
serialize("F.dat", F)
serialize("Zatoms.dat", Zatoms)
