module MyACE1pack

# load and reexport JuLIP, ACE1
# also ArgParse for the command line script
using Reexport 
@reexport using JuLIP
@reexport using ACE1
@reexport using ACE1x
@reexport using ACEfit
@reexport using ArgParse
export JuLIP, ACE1, ACEfit, Argparse

# Convenience Layer 

#include("atoms_data.jl")






include("fit.jl")

include("data.jl")

include("basis.jl")

include("model.jl")

include("solver.jl")

include("regularizer.jl")

include("read_params.jl")

include("export.jl")

include("analysis.jl")

# a little hack to load ACE1pack artifacts from anywhere? 
using LazyArtifacts
artifact(str) = (@artifact_str str)

end
