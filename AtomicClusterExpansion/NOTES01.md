  Presence of `Project.toml` in current directory will prevent
  use of `LOAD_PATH` to load module defined in a `.jl` file.

# Installation

ACESUITE environment is installed as shared environment.
You need to add ACEsuite registry first.

```
] activate @ACESUITE
```
or
```
Pkg.activate("ACESUITE, shared=true)
```

# Initializing basis functions

```julia
import ACE1x
basis = ACE1x.ace_basis(
    elements = [:Si,],
    order = 3,   
    totaldegree = 10,
    rcut = 5.0
)
```
The type of `basis` is `JuLIP.MLIPs.IPSuperBasis{JuLIP.MLIPs.IPBasis}`.

IPBasis seems to contain two basic elements:
- PolyPairBasis
- RPIBasis

