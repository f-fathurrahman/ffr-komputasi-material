

export acemodel

import MyJuLIP: energy, forces, virial, cutoff
import MyACE1.Utils: get_maxn
using LinearAlgebra: Diagonal 

_mean(x) = sum(x) / length(x)


"""
`struct MyACE1Model` : this specifies an affine `MyACE1.jl` model via a basis, the 
parameters and a reference potential `Vref`. The generic constructor is 
```
MyACE1Model(basis, params, Vref) 
``` 
where `basis` is typically a `MyACE1.RPIBasis` or a `MyJuLIP.MLIPs.IPSuperBasis` 
object, `params` is a vector of model parameters and `Vref` a reference 
potential, typically a `MyJuLIP.OneBody`. Setting new parameters is done 
via `MyACE1pack.set_params!`. 
A convenience constructor that exposes the most commonly used options to 
construct MyACE1 models is provided by `MyACE1pack.acemodel`. 
"""
mutable struct MyACE1Model 
   basis 
   params 
   Vref  
   potential
   meta::Dict
end

# TODO: 
# - some of this replicates functionality from MyACE1.jl and we should 
#   consider having it in just one of the two places. 
# - better defaults for transform
# - more documentation, especially if this becomes the standard interface  

"""
`function acemodel` : convenience constructor for `MyACE1Model` objects.
It takes only keyword arguments, are passed to `acebasis`. Please see the 
documentation of `acebasis` for the details. 

In addition to the `acebasis` arguments, `acemodel` also accepts the following
additional arguments: 
* `Eref` : reference energies for the species
* `Vref` : reference potential
If neither are provided then no reference potential is used. If both are provided 
then `Eref` is ignored. If only `Eref` is provided then a corresponding 
`Vref` one-body potential is constructed from the reference energies 
provided. 
"""
function acemodel(;  kwargs...)
   Eref = get(kwargs, :Eref, nothing)
   Vref = get(kwargs, :Vref, nothing)
   
   println("Eref = ", Eref)
   println("Vref = ", Vref)

   # construct the basis 
   basis = ace_basis(; kwargs...)

   if Vref == nothing && Eref != nothing 
      Vref = MyJuLIP.OneBody(Eref...)
   end

   # construct a model without parameters and without an evaluator 
   model = MyACE1Model(basis, nothing, Vref, nothing, Dict())

   # set some random parameters -> this will also generate the evaluator 
   params = randn(length(basis))
   params = params ./ (1:length(basis)).^4
   _set_params!(model, params)
   return model 
end


_sumip(pot1, pot2) = 
      MyJuLIP.MLIPs.SumIP([pot1, pot2])

_sumip(pot1::MyJuLIP.MLIPs.SumIP, pot2) = 
      MyJuLIP.MLIPs.SumIP([pot1.components..., pot2])
                                      
function _set_params!(model, params)
   model.params = params
   model.potential = MyJuLIP.MLIPs.combine(model.basis, model.params)
   if model.Vref != nothing 
      model.potential = _sumip(model.potential, model.Vref)
   end
   return model 
end



smoothness_prior(model::MyACE1Model; kwargs...) = 
      smoothness_prior(model.basis; kwargs...)
