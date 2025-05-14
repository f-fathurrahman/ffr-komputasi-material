```
python load_qm9.py dataset_qm9
```

Train QM9 energy (U0) prediction:

```
python train_energy_force.py dataset_qm9/qm9.db ./modeldir ./split50k.npz --ntrain 50000 --nval 10000 --fit_energy --atomref dataset_qm9/atomref.npz
```

