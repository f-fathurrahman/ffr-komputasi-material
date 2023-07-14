import JuLIP
import ACE1x
import Random
import ACEfit
import ACE1pack
using LinearAlgebra: Diagonal
import PrettyTables

include("MyACE1Model.jl")


# ### Step 1: specify the ACE Model
#
# The parameters have the following meaning: 
# * `elements`: list of chemical species, symbols 
# * `order` : correlation order
# * `totaldegree`: maximum total polynomial degree used for the basis 
# * `rcut` : cutoff radius (optional, defaults are provided)

model = my_acemodel(
    elements = [:Si,],
    order = 3,   
    totaldegree = 10,
    rcut = 5.0)
@show length(model.basis);