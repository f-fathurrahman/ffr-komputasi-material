using MyACEpotentials

data = read_extxyz("datasets/TiAl_tutorial.xyz")

# The next step is to generate a model
# order = 3 : We take 3-correlation, i.e. a 4-body potential,
# totaldegree = 6 : a very low polynomial degree just for testing
# rcut = 5.5 : this is a typical cutoff radius for metals
model = acemodel(
    elements = [:Ti, :Al],
    order = 3, totaldegree = 6,
    rcut = 5.5,
    Eref = [:Ti => -1586.0195, :Al => -105.5954]
)
println("Number of basis = ", length(model.basis))

# Specify weights used in loss function
weights = Dict(
    "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 1.0),
    "TiAl_T5000" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0)
)

# To estimate the parameters we still need to choose a solver for the least squares system.
# In this tutorial we use the LSQR algorithm, for no specific reason at all.
# Many other solvers are available, and can be explored by looking at the
# documentation of ACEfit.jl.
solver = MyACEfit.LSQR(damp=1e-2, atol=1e-6);

# ACE1.jl has a heuristic smoothness prior built in which assigns
# to each basis function Bi a scaling parameter si that estimates
# how "rough" that basis function is. The following line generates
# a regularizer (prior) with si^q on the diagonal, thus penalizing
# rougher basis functions and enforcing a smoother fitted potential.
P = smoothness_prior(model; p=4)    #  (p = 4 is in fact the default)

# Split dataset
train_data = data[1:5:end]
test_data = data[2:10:end]

acefit!(model, train_data; solver=solver, prior=P)

@info("Train Error Table")
MyACEpotentials.linear_errors(train_data, model)

@info("Test Error Table")
MyACEpotentials.linear_errors(test_data, model; weights=weights)

