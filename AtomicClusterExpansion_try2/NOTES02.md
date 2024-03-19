# Investigating ACE basis

The function: `ACE1x.ace_basis`.

Return type is quite long: `JuLIP.MLIPs.IPSuperBasis{JuLIP.MLIPs.IPBasis}`
This type only has one field: `BB` which is of type `Vector{IPBasis}`

The basis itself composed of the following types:

```julia
# julia> typeof(basis.BB[1])
PolyPairBasis{MyACE1.OrthPolys.TransformedPolys{Float64, MyACE1.Transforms.Agnesi2Transform{Float64, Int64}, MyACE1.OrthPolys.OrthPolyBasis{Float64}, PolyEnvelope{Float64}}, 2}
```
This is initialized by call to `ACE1x._pair_basis(kwargs)`
`PolyPairBasis` is defined in `ACE1.PairPotentials`.


```julia
# julia> typeof(basis.BB[2])
RPIBasis{Float64, BasicPSH1pBasis{Float64, 2, MyACE1.OrthPolys.TransformedPolys{Float64, MyACE1.Transforms.MultiTransform{2, MyACE1.Transforms.AffineT{Float64, MyACE1.Transforms.Agnesi2Transform{Float64, Int64}}}, MyACE1.OrthPolys.OrthPolyBasis{Float64}, MyACE1.OrthPolys.OneEnvelope}}, 2, MyACE1.DAGEvaluator}
```
This is initialized by call to `ACE1x.mb_ace_basis(kwargs)`
`RPIBasis` is defined in `ACE1`.


## Pair Basis

Type of pair basis is `MyACE1.PairPotentials.PolyPairBasis` which full names are quite
long (it includes type parameters).
The fields of PolyPairBasis:
```julia
# typeof(pairB) |> fieldnames
(:zlist, :J, :bidx0)
```

`zlist` contains species information (their atomic numbers)

`bidx0` is a `SMatrix` of size (Nspecies,Nspecies) containing indices (of what???).

`J` is the actual pair basis, which is a `SMatrix` of size (Nspecies,Nspecies).

The element type of `J` is `ACE1.OrthPolys.TransformedPolys`.

```julia
# pairB.J |> eltype
MyACE1.OrthPolys.TransformedPolys{Float64, MyACE1.Transforms.Agnesi2Transform{Float64, Int64}, MyACE1.OrthPolys.OrthPolyBasis{Float64}, MyACE1.OrthPolys.PolyEnvelope{Float64}}
```

The fields of `TransformedPolys` are

```julia
# julia> pairB.J[1,1] |> typeof |> fieldnames
(:J, :trans, :rl, :ru, :envelope)
```

There is `J` again.
```julia
# julia> pairB.J[1,1].J |> typeof
MyACE1.OrthPolys.OrthPolyBasis{Float64}

# julia> pairB.J[1,1].trans |> typeof
MyACE1.Transforms.Agnesi2Transform{Float64, Int64}

# julia> pairB.J[1,1].rl |> typeof
Float64

# julia> pairB.J[1,1].ru |> typeof
Float64

# julia> pairB.J[1,1].envelope |> typeof
MyACE1.OrthPolys.PolyEnvelope{Float64}
```

Fields of `OrthPolyBasis` are simple.

`evaluate!` is imported from JuLIP.


## MB basis (?)





# Some questions

How to use spherical harmonics basis:
  module ACE1.SphericalHarmonics

How to use radial basis

How to use 1p basis
how to check permutation invariance
how to check rotation invariance

BasicPSH1pBasis

Pure2b ? What is pure here means?