`qeManager` adalah pustaka Python yang dapat digunakan untuk
membuat input untuk Quantum ESPRESSO dan beberapa paket komputasi
*ab-initio* yang lain.

```python
import sys
sys.path.append('/home/efefer/WORKS/my_github_repos/')
```

## Manipulating PWSCF input directly

Create an input for NH3 in box

First we import several modules:

```python
from ase import Atoms
from ase.units import Bohr
import numpy as np
from qeManager.pwscf import *
from ase.build import molecule
```

We build NH3 molecule, using ASE

```python
atoms = molecule("NH3")
atoms.set_pbc([True,True,True])
cell = np.array([16.0,16.0,16.0])*Bohr
atoms.set_cell(cell)
atoms.center()
```
Now we can start building an input for PWSCF.

We start by `CONTROL` namelist:
```python
ctrl_NL = ControlNameList() # using default parameters
```
The above constructor will initialize `CONTROL` namelist with default
parameters. Several parameters still need to be set manually.

For example we set `pseudo_dir` manually.

```python
ctrl_NL.pseudo_dir = "/home/efefer/pseudo"
```

Write all parameters, including defaults to stdout

```python
ctrl_NL.write_all()
```

Also do the same for `SYSTEM` and `ELECTRONS`

```python
sys_NL = SystemNameList(atoms)
sys_NL.write_all()

elec_NL = ElectronsNameList()
elec_NL.write_all()
```

```python
# pseudopotentials list
pspFiles = ["N_ONCV_PBE-1.0.upf", "H_ONCV_PBE-1.0.upf"]
# write ATOMIC_SPECIES card
write_atomic_species(atoms, pspFiles=pspFiles)
# write ATOMIC_POSITIONS
write_atomic_positions(atoms)
# write CELL_PARAMETERS
write_cell(atoms)
```