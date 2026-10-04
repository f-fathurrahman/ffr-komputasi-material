# Machine Learning XC (ML XC)

This module provides the modern ML XC pipeline used by `atom`, including
model wrappers, preprocessing metadata, and SCF integration.

## Key components

- `MLXCCalculator`: wraps a trained model together with preprocessing metadata
  (scalers, symlog settings, target functional/component/mode).
- `TorchXCModel`: torch backend that handles training, evaluation, and I/O.
- `VxcDataLoader` / `ExcDataLoader`: data containers produced by the data
  pipeline and consumed by `TorchXCModel` and `MLXCCalculator`.
- `nn_models`: optional torch model definitions (e.g., `MLP`,
  `ChannelEmbeddingMLP`, `ChannelEmbeddingResNet`).

## Supported model kinds

- **`potential`**: predicts XC potentials (e.g., `v_xc`, `v_x`, `v_c`).
  Fully supported end-to-end.
- **`energy`**: placeholders exist, but energy training/evaluation and
  `MLXCCalculator.compute_energy()` are not implemented yet.

## Valid feature sets

These are validated by `data_loading.py`:

- **Potential models**: `["rho", "grad_rho", "lap_rho", "hartree", "lda_xc"]`
- **Energy density models**: `["rho", "grad_rho", "lap_rho"]`

## Metadata and file layout

`MLXCCalculator.save()` produces the following files:

- `<name>_model.pth`
- `<name>_config.json`
- `<name>_scalers.pkl`
- `<name>_metadata.json`

## Typical workflow

### 1) Train a potential model

```python
from atom.xc.ml_xc import TorchXCModel, MLXCCalculator
from atom.xc.ml_xc.nn_models import MLP

# train_loader / val_loader are VxcDataLoader instances
model = TorchXCModel(
    model_kind="potential",
    features_list=["rho", "grad_rho", "lap_rho", "hartree", "lda_xc"],
    model_cls=MLP,
    model_init_kwargs={"input_dim": 5, "output_dim": 1},
    model_name="pbe_delta",
    model_dir="./models",
    device="cpu",
)

model.train(train_loader, val_loader=val_loader, epochs=100, batch_size=256)

mlxc = MLXCCalculator.from_dataloader(model=model, data_loader=train_loader)
mlxc.save("./models", overwrite=True)
```

### 2) Load a saved calculator

```python
from atom.xc.ml_xc import MLXCCalculator

mlxc = MLXCCalculator.load("./models", "pbe_delta", device="cpu")
mlxc.print_info()
```

### 3) SCF integration (delta learning)

When `ml_xc_calculator` is provided, the solver expects
`ml_xc_calculator.target_functional == xc_functional`, then switches to
`ml_xc_calculator.reference_functional` internally and applies the ML delta.

```python
from atom.solver import AtomicDFTSolver

solver = AtomicDFTSolver(
    atomic_number=8,
    xc_functional="GGA_PBE",
    ml_xc_calculator=mlxc,
)
```

## Notes

- `target_component` is typically one of: `v_xc`, `v_x`, `v_c`, `v_x_v_c`.
- `target_mode` is either `absolute` or `delta`.
- If `model_kind == "potential"`, ML energy correction is not available.

