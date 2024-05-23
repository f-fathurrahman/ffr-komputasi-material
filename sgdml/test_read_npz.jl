using NPZ

F = npzread("DATASET/ethanol_dft.npz", ["F"])["F"]
R = npzread("DATASET/ethanol_dft.npz", ["R"])["R"]
E = npzread("DATASET/ethanol_dft.npz", ["E"])["E"]
# Convert to Int64
z = Int64.( npzread("DATASET/ethanol_dft.npz", ["z"])["z"] )



