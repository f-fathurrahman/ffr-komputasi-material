using Serialization: deserialize

R = deserialize("R.dat")
E = deserialize("E.dat")
F = deserialize("F.dat")
Zatoms = deserialize("Zatoms.dat")

Ndata = length(E)
@assert Ndata == length(R)
@assert Ndata == length(F)

Natoms = length(Zatoms)
@assert Natoms == size(R[1],2)
@assert Natoms == size(F[1],2)


dim_i = 3*Natoms
desc_dim = (Natoms * (Natoms - 1)) / 2 |> Int64


#tril_indices = np.tril_indices(n_atoms, k=-1)
