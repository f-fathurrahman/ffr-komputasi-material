# Investigating ACE basis

The function: `ACE1x.ace_basis`.

Return type is quite long: `JuLIP.MLIPs.IPSuperBasis{JuLIP.MLIPs.IPBasis}`
This type only has one field: `BB` which is of type `Vector{IPBasis}`

The basis itself composed of the following types:

```julia
# julia> typeof(basis.BB[1])
PolyPairBasis{...}
```
This is initialized by call to `ACE1x._pair_basis(kwargs)`
`PolyPairBasis` is defined in `ACE1.PairPotentials`.


```julia
# julia> typeof(basis.BB[2])
RPIBasis{...}
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
MyACE1.OrthPolys.TransformedPolys{...}
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

To construct MB basis, we first need to create radial basis `rbasis`
```julia
# julia> typeof(rbasis)
MyACE1.OrthPolys.TransformedPolys{...}
```

Then, using Pure2b we can create rpibasis. Our kwargs lead to Pure2b.

Fields of RPIBasis:
```
julia> rpibasis |> show_fields

 Type of variable: MyACE1.RPI.RPIBasis{....}

  fieldname = pibasis, type = MyACE1.PIBasis{...}

  fieldname = A2Bmaps, type = Tuple{SparseArrays.SparseMatrixCSC{Float64, Int64}, SparseArrays.SparseMatrixCSC{Float64, Int64}}

  fieldname = Bz0inds, type = Tuple{UnitRange{Int64}, UnitRange{Int64}}
```

A2Bmaps[1] and A2Bmaps[2] are the same ?

It seems that the actual basis is in `pibasis`

```
julia> rpibasis.pibasis |> show_fields

 Type of variable: MyACE1.PIBasis{...}

  fieldname = basis1p, type = MyACE1.RPI.BasicPSH1pBasis{...}

  fieldname = zlist, type = MyJuLIP.Potentials.SZList{2}

  fieldname = inner, type = Tuple{MyACE1.InnerPIBasis, MyACE1.InnerPIBasis}

  fieldname = evaluator, type = MyACE1.DAGEvaluator
```



Fields of `rpibasis.pibasis.basis1p`
```
julia> rpibasis.pibasis.basis1p |> show_fields

 Type of variable: MyACE1.RPI.BasicPSH1pBasis{...}

  fieldname = J, type = MyACE1.OrthPolys.TransformedPolys{...}

  fieldname = SH, type = MyACE1.SphericalHarmonics.SHBasis{Float64}

  fieldname = zlist, type = MyJuLIP.Potentials.SZList{2}

  fieldname = spec, type = Vector{MyACE1.RPI.PSH1pBasisFcn}

  fieldname = Aindices, type = Matrix{UnitRange{Int64}}
```

Is rpibasis.pibasis.basis1p.J the same as `rbasis`?

Fields of `rpibasis.pibasis.basis1p.J`:
```
julia> rpibasis.pibasis.basis1p.J |> show_fields

 Type of variable: MyACE1.OrthPolys.TransformedPolys{...}

  fieldname = J, type = MyACE1.OrthPolys.OrthPolyBasis{Float64}

  fieldname = trans, type = MyACE1.Transforms.MultiTransform{...}

  fieldname = rl, type = Float64

  fieldname = ru, type = Float64

  fieldname = envelope, type = MyACE1.OrthPolys.OneEnvelope
```

Now, we get `OrthPolyBasis`:
```
julia> rpibasis.pibasis.basis1p.J.J |> show_fields

 Type of variable: MyACE1.OrthPolys.OrthPolyBasis{Float64}

  fieldname = pl, type = Int64

  fieldname = tl, type = Float64

  fieldname = pr, type = Int64

  fieldname = tr, type = Float64

  fieldname = A, type = Vector{Float64}

  fieldname = B, type = Vector{Float64}

  fieldname = C, type = Vector{Float64}

  fieldname = tdf, type = Vector{Float64}

  fieldname = ww, type = Vector{Float64}
```


```
julia> rpibasis.pibasis.basis1p.J.trans |> show_fields

 Type of variable: MyACE1.Transforms.MultiTransform{...}

  fieldname = zlist, type = MyJuLIP.Potentials.SZList{2}

  fieldname = transforms, type = StaticArraysCore.SMatrix{...}
```

How more informations are added to `OrthPolyBasis` until we get `RPIBasis` ?


# Some questions

How to use spherical harmonics basis:
  module ACE1.SphericalHarmonics

How to use radial basis

How to use 1p basis
how to check permutation invariance
how to check rotation invariance

BasicPSH1pBasis

Pure2b ? What is pure here means?