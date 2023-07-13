  Presence of `Project.toml` in current directory will prevent
  use of `LOAD_PATH` to load module defined in a `.jl` file.

ACESUITE environment is installed as shared environment.
You need to add ACEsuite registry first.

```
] activate @ACESUITE
```
or
```
Pkg.activate("ACESUITE, shared=true)
```


