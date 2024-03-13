module MyACEbase

using Reexport, StaticArrays

# @def macro generator
include("def.jl")

# evaluation interface functions to be shared across ACEsuit packages. the 
# actual interfaces for these can vary across packages and is not 
# restricted in any way. 
function evaluate end
function evaluate_d end
function evaluate_dd end
function evaluate_ed end
function evaluate_ed2 end
function evaluate! end
function evaluate_d! end
function evaluate_dd! end
function evaluate_ed! end
function evaluate_ed2! end

# these could be useful functions to share across ACEsuit packages as well since 
# there are many different ways how one can allocate and release memory. 
# For now, we will make them owned by ObjectPools.jl so that that package 
# doesn't have to depend on MyACEbase. 
# function acquire! end 
# function release! end 

# This is a helpful little utility that we don't know where else to put.
"""
a simple utility function to check whether two objects are equal
"""
allfieldsequal(x1, x2) =
      all( getfield(x1, sym) == getfield(x2, sym)
           for sym in union(fieldnames(typeof(x1)), fieldnames(typeof(x2))) )

@deprecate _allfieldsequal(args...) allfieldsequal(args...)

# This creates a sub-module FIO
# FIO via JSON and YAML functionality is collected here. We don't trust 
# any of the standard packages and write our own manual object <-> JSON
# routines. 
include("fio.jl")
@reexport using MyACEbase.FIO

# This creates a sub-module Testing which several nice testing routines
# that can be shared across packages.  
include("testing.jl")

end
