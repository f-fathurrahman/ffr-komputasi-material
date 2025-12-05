using MyACEpotentials

data = read_extxyz("datasets/TiAl_tutorial.xyz")

#=
The next step is to generate a basis set:
- order = 3 : We take 3-correlation, i.e. a 4-body potential,
- totaldegree = 6 : a very low polynomial degree just for testing
- rcut = 5.5 : this is a typical cutoff radius, there is also a good
  default which is a bit higher

These three are the most important approximation parameters to explore
when trying to improve the fit-accuracy. In addition there is
- The parameter r0 is just a scaling parameter and the fits should not be very
  sensitive to its choice. A rough estimate for the nearest-neighbour distance
  is usually ok. (NB: if you change to a non-trivial distance transform,
  then the parameter r0 may become important.)
=#

r0 = 2.88
basis = MyACE1x.ace_basis(
    elements = [:Ti, :Al],
    order = 3,
    totaldegree = 6,
    rcut = 5.5,
    r0 = r0)
@show length(basis);
# r0: specific to basis only?
# In acemodel function what is the default value of r0?
# We don't specify Eref or Vref for ace_basis, this is specified in OneBody

# Vref specifies a reference potential, which is subtracted from the training data and
# the ACE parameters are then estimated from the difference. This reference potential
# will in the end be added to the ACE model. Here we use a one-body potential i.e.
# a reference atom energy for each individual species. Usage of a one-body reference
# potential generally results in very slightly reduced fit accuracy but significantly
# improved 2-body potentials with a realistic dimer shape.

Vref = OneBody(:Ti => -1586.0195, :Al => -105.5954);

# For loss function
weights = Dict(
        "FLD_TiAl" => Dict("E" => 60.0, "F" => 1.0 , "V" => 1.0 ),
        "TiAl_T5000" => Dict("E" => 5.0, "F" => 1.0 , "V" => 1.0 ));
# What is the rationale for these values?


datakeys = (energy_key = "energy", force_key = "force", virial_key = "virial")

train_data = data[1:5:end]
# Construct training data as AtomsData
train = [
    MyACEpotentials.AtomsData(t; weights=weights, v_ref=Vref, datakeys...) for t in train_data
]

# Assemble matrices
A, Y, W = MyACEfit.assemble(train, basis);


# ACE1.jl has a heuristic smoothness prior built in which assigns to each basis
# function Bi a scaling parameter si that estimates how "rough" that basis function is.
# The following line generates a regularizer (prior) with si^q on the diagonal,
# thus penalizing rougher basis functions and enforcing a smoother fitted potential.
P = smoothness_prior(basis; p = 4)

# Once all the solver parameters have been determined, we use ACEfit
# to estimate the parameters. This routine will return the fitted interatomic
# potential IP as well as the a dictionary lsqfit with some information about the fitting process.

# Create solver
solver = MyACEfit.LSQR(damp = 1e-2, atol = 1e-6, P = P)

# Solve the linear problem
results = MyACEfit.solve(solver, W .* A, W .* Y)

# Construct a JuLIP potential from basis and its coefficients
pot_1 = MyJuLIP.MLIPs.SumIP(Vref, MyJuLIP.MLIPs.combine(basis, results["C"]));



# The advantage of working with the ACE basis rather than the ACE model interface
# is that we can now make some changes to the fitting parameters and refit.
# For example, we might want different weights, change the smoothness prior,
# and switch to a RRQR solver.

weights["FLD_TiAl"]["E"] = 20.0
W = MyACEfit.assemble_weights(train)
# ??? `weights` also changes in `train` ?

new_prior = smoothness_prior(basis; p = 2)
solver = MyACEfit.RRQR(; rtol=1e-8, P=new_prior)

results = MyACEfit.solve(solver, W .* A, W .* Y)
# Construct another JuLIP basis
pot_2 = MyJuLIP.MLIPs.SumIP(Vref, MyJuLIP.MLIPs.combine(basis, results["C"]));