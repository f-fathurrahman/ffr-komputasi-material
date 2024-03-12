# MyJuLIP.jl master file.

module MyJuLIP

using Reexport
@reexport using NeighbourLists

const _usethreads = Ref(true)
function usethreads!(tf::Bool)
   MyJuLIP._usethreads[] = tf
end
nthreads() = MyJuLIP._usethreads[] ? Threads.nthreads() : 1

# quickly switch between Matrices and Vectors of SVectors, etc
include("arrayconversions.jl")

import ACEbase: FIO
@reexport using MyJuLIP.FIO

# define types and abstractions of generic functions
include("abstractions.jl")

include("chemistry.jl")
@reexport using MyJuLIP.Chemistry

# the main atoms type
include("atoms.jl")
include("dofmanagement.jl")

# how to build some simple domains
include("build.jl")
@reexport using MyJuLIP.Build

# a few auxiliary routines
include("utils.jl")
@reexport using MyJuLIP.Utils

# interatomic potentials prototypes and some example implementations
include("Potentials.jl")
@reexport using MyJuLIP.Potentials
# and we want to import some more functions from `Potentials` which are really
# central to MyJuLIP, so that they can be extended using just `import MyJuLIP: ...`
import MyJuLIP.Potentials: numz, z2i, i2z

# basic preconditioning capabilities
include("preconditioners.jl")
@reexport using MyJuLIP.Preconditioners

# some solvers
include("Solve.jl")
@reexport using MyJuLIP.Solve

# experimental features
include("Experimental.jl")
@reexport using MyJuLIP.Experimental

# the following are some sub-modules that are primarily used
# to create further abstractions to be shared across several
# modules in the MyJuLIP-verse.
include("mlips.jl")


# codes to facilitate testing
include("Testing.jl")
include("datadeps.jl")

end # module
