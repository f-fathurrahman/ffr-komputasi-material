module MyACEpotentials

using Reexport 
@reexport using MyJuLIP
@reexport using MyACE1
@reexport using MyACE1x
@reexport using MyACEfit
@reexport using MyACEmd

include("atoms_data.jl")
include("model.jl")
include("export.jl")
#include("example_data.jl")
include("descriptor.jl")
include("atoms_base.jl")
include("io.jl")

include("analysis/potential_analysis.jl")
include("analysis/dataset_analysis.jl")

include("outdated/fit.jl")
include("outdated/data.jl")
include("outdated/basis.jl")
include("outdated/solver.jl")
include("outdated/regularizer.jl")
include("outdated/read_params.jl")

end
