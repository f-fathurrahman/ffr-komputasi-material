# Post-processing Utilities for Atomic DFT

This module provides utilities for atomic DFT workflows, including data generation, loading, processing, evaluation, and visualization. The primary entry point is `AtomicDataManager`, which coordinates the generator, loader, and processor.

## Usage

### Import the utilities

```python
from atom.data import (
    AtomicDataManager,
    DataGenerator,
    DataLoader,
    DataProcessor,
    DataVisualizer,
    AtomicDataset,
    VxcDataLoader,
    SingleConfigurationData,
)
```

### Example: Generate data

```python
manager = AtomicDataManager(
    data_root="data",
    scf_xc_functional="pbe0",
    forward_pass_xc_functionals=["gga_pbe"],
)

# This will prompt before overwriting existing data.
manager.generate_data(
    atomic_number_list=[1],  # Hydrogen
    domain_size=20.0,
    finite_elements_number=15,
    polynomial_order=20,
    quadrature_point_number=43,
    oep_basis_number=5,
    mesh_type="polynomial",
    mesh_concentration=2.0,
    mesh_spacing=0.1,
    start_configuration_index=1,  # Starting index for configuration folders (default: 1)
    overwrite=False,  # If True, skip confirmation prompts (default: False)
)
```

### Example: Load dataset and prepare a dataloader

```python
dataset = manager.load_data(
    atomic_number_list=None,  # Load all available atoms
    use_radius_cutoff=True,
    smooth_vxc=False,
    include_intermediate=False,
    print_summary=True,
)

vxc_loader = dataset.prepare_potential_dataloader(
    target_functional="gga_pbe",
    target_component="v_xc",
    reference_functional="pbe0",
)

X = vxc_loader.X
y = vxc_loader.y
weights = vxc_loader.weights_for_potential
atomic_numbers = vxc_loader.atomic_numbers_per_sample
```


### Example: Visualize results

```python
# Plot predictions vs true values
fig, ax = DataVisualizer.plot_predictions_vs_true(
    y_true=y_test,
    y_pred=y_pred,
    model_name="MyModel",
    target_name="Potential",
    atomic_numbers=atomic_numbers_test,  # Color points by atomic number
    save_path="predictions.png",
)

# Plot error vs radius
fig, ax = DataVisualizer.plot_difference_vs_radius(
    radius=radius_test,
    y_true=y_test,
    y_pred=y_pred,
    model_name="MyModel",
    atomic_numbers=atomic_numbers_test,
    save_path="error_vs_radius.png",
)
```

## Module Structure

- **data_generation.py**: Dataset generation utilities (`DataGenerator`)
- **data_loading.py**: Data loading utilities (`DataLoader`, `SingleConfigurationData`)
- **data_manager.py**: Dataset management and preparation (`AtomicDataset`, `AtomicDataManager`, `VxcDataLoader`)
- **data_processing.py**: Data processing utilities (`DataProcessor`)
- **data_visualization.py**: Visualization utilities (`DataVisualizer`)
- **evaluation.py**: Model evaluation utilities (`evaluate_model`, `print_metrics`)

## Notes

- Model evaluation utilities require `scikit-learn`.
- Visualization utilities require `matplotlib`.
- The module can be used both in Python scripts and Jupyter notebooks.

