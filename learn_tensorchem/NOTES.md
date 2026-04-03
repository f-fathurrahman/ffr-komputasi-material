Build extension (calling Fortran subroutines)
```
f2py3 -c maxvol.f90 -m maxvol -lopenblas
```
It is important to also specify link to BLAS.

