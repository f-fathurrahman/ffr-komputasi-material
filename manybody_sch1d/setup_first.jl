using Pkg
Pkg.activate("SCH1D", shared=true)

using Revise, Infiltrator
using SparseArrays, LinearAlgebra
using Arpack: eigs
using Combinatorics: permutations