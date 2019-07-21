I just noticed that PWSCF (as distributed in Quantum Espresson version 6.4)
has a new diagonalization method: PPCG (projected
preconditioned conjugate gradient). However, it not documented yet
in the HTML input documentation.

Now, PWSCF have three methods for iterative diagonalization of Hamiltonian:

- Davidson method (default)
- Sequential CG
- PPCG

In this post, I want make a comparison between these methods (especially
between Davidson and PPCG).

I consider two solid-state systems for this: silicon fcc and
LiNiO<sub>2</sub>.

**Si fcc**

First, let's consider silicon fcc. Here is the input file:

```
&CONTROL
  calculation = 'scf'
  pseudo_dir = '/home/efefer/pseudo'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 2
  celldm(1) = 10.20
  nat = 2
  ntyp = 1
  ecutwfc = 18.0
/

&ELECTRONS
  diagonalization = 'david'
  mixing_beta = 0.7
  conv_thr =  1.0d-8
/

ATOMIC_SPECIES
Si  28.086  Si.pz-vbc.UPF

ATOMIC_POSITIONS crystal
Si 0.00 0.00 0.00
Si 0.25 0.25 0.25

K_POINTS automatic
  8 8 8 0 0 0
```

To set the diagonalization to Davidson method, we use the line:
```
  diagonalization = 'david'
```

For CG method:
```
  diagonalization = 'cg'
```

And for PPCG method:
```
  diagonalization = 'ppcg'
```

Note that I have used tighter convergence criteria `ethr` than the default value.
Because this system is very small, I will only use 1 processor for the test.

The first thing that I want to compare is the memory requirement:
```
LOG_cg   :     Estimated max dynamical RAM per process >       2.53 MB
LOG_david:     Estimated max dynamical RAM per process >       2.63 MB
LOG_ppcg :     Estimated max dynamical RAM per process >       2.53 MB
```
We can see that PPCG method requires less memory than Davidson method. Its memory requirement
is about the same as CG method. This is an advantage of CG and PPCG as compared to
Davidson method.

Now let's see the time required to finish the calculation:
```
LOG_cg   :     PWSCF        :      0.27s CPU      0.32s WALL
LOG_david:     PWSCF        :      0.27s CPU      0.28s WALL
LOG_ppcg :     PWSCF        :      0.33s CPU      0.37s WALL
```
Here we see that the job PPCG requires the longest time to finish
as compare to the job that used CG or Davidson method.

Here we see that the job took less than 1 second to finish for the
three methods that we considered here.
This is not surprising since the system that we test, the silicon fcc, is
very small. This means that the timing should not be taken seriously and
the timing result might vary slightly between experiments (measurement).
Also, we also should take note that the overall time to finish the calculation
also depend on the other aspects of SCF algorithm such as density mixing and
not from the effect of diagonalization method only.
However, in the comparison test that I made, I only varied the
the value of `diagonalization` in the input file, so we may assume that
the difference comes from the diagonalization method.

Let's check how many SCF iterations it took to reach convergence:
```
LOG_cg   :     convergence has been achieved in   6 iterations
LOG_david:     convergence has been achieved in   6 iterations
LOG_ppcg :     convergence has been achieved in   5 iterations
```

We see that SCF+PPCG only requires 5 iterations to complete as opposed
to SCF+CG and SCF+Davidson which require 6 iterations. Even though so,
PPCG requires more time to finish the job (for this case).
From this result, I guess that the overhead time (such as initialization, etc)
for PPCG is larger as compared to CG and Davidson. This is understandable
since it is only recently developed.

Let's not forget to check the final total energy result:
```
LOG_cg   :!    total energy              =     -15.84445125 Ry
LOG_david:!    total energy              =     -15.84445125 Ry
LOG_ppcg :!    total energy              =     -15.84445125 Ry
```
The diagonalization method should not affect the converged total energy
result.

**LiNiO<sub>2</sub> using PAW dataset**

Now let's consider a larger method (but still manageable on a laptop).
Here an example input file:

```
&CONTROL
  calculation = 'scf'
  pseudo_dir = '/home/efefer/pseudo'
  outdir = './tmp'
/

&SYSTEM
  ibrav = 0
  A = 5.07202
  nat = 4
  ntyp = 3
  ecutwfc = 40.0
  ecutrho = 180.0
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.002
/

&ELECTRONS
  diagonalization = 'david'
  electron_maxstep = 150
  mixing_beta = 0.1
/

ATOMIC_SPECIES
  Ni   58.69340  Ni.pbe-n-kjpaw_psl.0.1.UPF
   O   15.99900  O.pbe-n-kjpaw_psl.0.1.UPF
  Li    6.96750  Li.pbe-s-kjpaw_psl.0.2.1.UPF

ATOMIC_POSITIONS {crystal}
Ni   0.000000000000000   0.000000000000000   0.000000000000000 
 O   0.741937000000000   0.741937000000000   0.741937000000000 
 O   0.258063000000000   0.258063000000000   0.258063000000000 
Li   0.500000000000000   0.500000000000000   0.500000000000000 

K_POINTS automatic
5 5 5  0 0 0

CELL_PARAMETERS {alat}
  1.000000000000000   0.000000000000000   0.000000000000000 
  0.836301242074196   0.548270218510140   0.000000000000000 
  0.836301242074196   0.249697083586569   0.488110232379448 
```

Note that I have used PAW method for this system. The cutoff and kpoints parameters
are not optimized, I only used them for convenience.
The value of `ethr` is set to the default value (10<sup>-6</sup> Ry)

First, let's check the total energy result:
```
LOG_cg   :!    total energy              =    -284.13796750 Ry
LOG_david:!    total energy              =    -284.13796726 Ry
LOG_ppcg :!    total energy              =    -284.13796748 Ry
```
The values are agree within `ethr`.

Next, let's check the estimated RAM usage:
```
LOG_cg   :     Estimated max dynamical RAM per process >      61.96 MB
LOG_david:     Estimated max dynamical RAM per process >      61.96 MB
LOG_ppcg :     Estimated max dynamical RAM per process >      61.96 MB
```
Hmm, they are similar. This is actually outside my expectation as I expect
they should have difference (at least for CG and Davidson). I have no explanation
or guess other than that they are not estimated properly (?).

Next, let's check the total time required to finish the calculation:
```
LOG_cg   :     PWSCF        :   2m43.99s CPU   2m45.24s WALL
LOG_david:     PWSCF        :     54.25s CPU     54.82s WALL
LOG_ppcg :     PWSCF        :   2m 3.64s CPU   2m 4.19s WALL
```
CG is the slowest and Davidson is the fastest. PPCG is faster several seconds
as compared to CG.

Finally, let's check the total number SCF of iterations required
to reach convergence:
```
LOG_cg:     convergence has been achieved in  23 iterations
LOG_david:     convergence has been achieved in  12 iterations
LOG_ppcg:     convergence has been achieved in  12 iterations
```

SCF+Davidson and SCF+PPCG require the same number of iterations
to converges.
SCF+CG requires more iterations to converge however, if we compare the
time require for each SCF iterations:
```
In [1]: time_cg = 2*60 + 43.99;

In [2]: time_ppcg = 2*60 + 3.64;

In [3]: time_cg
Out[3]: 163.99

In [4]: time_ppcg
Out[4]: 123.64

In [5]: time_ppcg/12
Out[5]: 10.303333333333333

In [6]: time_cg/23
Out[6]: 7.130000000000001
```
Average time required per SCF iteration for CG is still faster than PPCG.
What I can conlude from this is that the resulting eigenpairs
of PPCG is more accurate than CG. Usually, we don't need the resulting eigenpairs
in each SCF iterations to be accurate, especially in early SCF iterations, as long
as they reach full accuracies in the last SCF iteration.
From my experience in [PWDFT.jl](https://github.com/f-fathurrahman/PWDFT.jl),
a good convergence in eigenpairs usually will result in faster SCF convergence
(less iterations to reach convergence).


**LiNiO<sub>2</sub> using ONCV pseudopotentials**

Because I still suspect the RAM requirement result from the previous
test, I decided to use ONCV pseudopotential,
which is a type of norm-conserving pseudopotential, instead of PAW.

```

```

