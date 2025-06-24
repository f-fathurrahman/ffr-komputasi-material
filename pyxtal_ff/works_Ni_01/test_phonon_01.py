import sys
sys.path.append('../')

from mini_mlip import MiniMLIP, MiniMLIPCalculator

train_data = "ALL_ATOMS.xyz"
path_model = "LOGDIR_mini_Ni_fcc_352/" # need trailing /

descriptor_dict = {
    "Rc": 4.0,
    "type": "SO3",
    "parameters": {
        "nmax": 4,
        "lmax": 3
    },
    "ncpu": 1,
}

model_dict = {
    "algorithm": "PR",
    "system" : ["Ni"],
    "path": path_model,
    'force_coefficient': 0.001,
}

ff = MiniMLIP(model=model_dict, descriptors=descriptor_dict)
ff.run(mode='predict', mliap=path_model+"PolyReg-checkpoint.pth")
calc = MiniMLIPCalculator(ff=ff)

from ase.build import bulk
from my_phonons import Phonons

# Setup crystal and EMT calculator
atoms = bulk("Ni", "fcc")

# Phonon calculator
N = 7
ph = Phonons(atoms, calc, supercell=(N, N, N), delta=0.05)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
ph.clean()

path = atoms.cell.bandpath("GXULGK", npoints=100)
bs = ph.get_band_structure(path)

dos = ph.get_dos(kpts=(20, 20, 20)).sample_grid(npts=100, width=1e-3)

# Plot the band structure and DOS:
import matplotlib.pyplot as plt
plt.style.use("dark_background")

fig = plt.figure(1, figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.35
bs.plot(ax=ax, emin=0.0, emax=emax)

dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color="grey",
                   edgecolor="k", lw=1)

dosax.set_ylim(0, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)

fig.savefig("TEMP_Ni_phonon.png", dpi=150)

