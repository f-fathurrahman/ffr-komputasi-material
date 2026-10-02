## Loading data

Load data, e.g. from VASP:
```julia
datadir = "./data_Si_vasp"
data = load_vasp(datadir)
```
Only need to specify `datadir`.


```sh
ls data_Si_vasp/
POSCAR  vasprun.xml
```

Some quantities in `data`:
```
julia> fieldnames(typeof(data))
(:lattice, :positions, :species, :kpoints, :weights, :ebands, :occupations, :fermi, :nelect, :magmom)
```