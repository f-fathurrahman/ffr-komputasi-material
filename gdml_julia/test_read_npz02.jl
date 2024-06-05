using NPZ

idxs_train = npzread("idxs_train_strat_sample.npz")["arr_0"] .+ 1
# need to offset with 1 because these indices are originally from Python

using Serialization
serialize("idxs_train.dat", idxs_train)
println("File serialized")