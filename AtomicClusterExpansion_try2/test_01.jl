include("setup_works.jl")

using MyACEpotentials: acemodel
import Random
using LinearAlgebra: norm, Diagonal

model = acemodel(elements=[:Si,], order=3, totaldegree=10, rcut=5.0)

