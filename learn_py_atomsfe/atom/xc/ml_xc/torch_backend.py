"""
Torch backend for XC models.
"""

from __future__ import annotations

import numpy as np
import json
import warnings
from pathlib import Path
from dataclasses import field, asdict
from typing import Any, Dict, Optional, Type, List, Literal, Union, Tuple


try:
    import torch
    from torch import nn
    TORCH_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    torch = None
    nn = None
    TORCH_AVAILABLE = False

try:
    import sklearn
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


from .base import XCBaseModel, ModelKind, ConfigExt, DataLoaderType, EvaluationMetrics, PotentialEvaluationMetrics, EnergyEvaluationMetrics
from ...data.data_manager import VxcDataLoader, ExcDataLoader, ComputationConfig




MODEL_NOT_NNMODULE_ERROR = \
    "Loaded model is not a torch.nn.Module, get {} instead."
DATA_LOADER_NOT_VXCDATALOADER_ERROR = \
    "parameter 'data_loader' must be a VxcDataLoader, get {} instead."
DATA_LOADER_NOT_EXCDATALOADER_ERROR = \
    "parameter 'data_loader' must be a ExcDataLoader, get {} instead."

# Training parameters error messages
FEATURE_NOT_IN_DATA_LOADER_ERROR = \
    "feature {} is not in data loader's features list {}. Check the features list of the data loader."
VXC_DATA_LOADER_CONTAINS_NO_SAMPLES_ERROR = \
    "VxcDataLoader contains no samples, please check the data loader to be non-empty."
EXC_DATA_LOADER_CONTAINS_NO_SAMPLES_ERROR = \
    "ExcDataLoader contains no samples, please check the data loader to be non-empty."
OPTIMIZER_NOT_VALID_ERROR = \
    "Invalid optimizer: {}. Use 'Adam' or 'SGD'."
CRITERION_NOT_VALID_ERROR = \
    "Invalid criterion: {}. Use 'L1Loss', 'MSELoss', 'L1norm', or 'L2norm'."
DATA_LOADER_TYPE_NOT_VALID_ERROR = \
    "Data loader type is not valid: {}. Use VxcDataLoader or ExcDataLoader."
WEIGHTS_NOT_SAME_AS_SAMPLES_ERROR = \
    "The number of weights must be the same as the number of samples, get {} weights and {} samples."
TORCH_NOT_AVAILABLE_FOR_TORCHXCMODEL_ERROR = \
    "TorchXCModel requires torch to be installed."
SKLEARN_NOT_AVAILABLE_FOR_MODEL_EVALUATION_ERROR = \
    "TorchXCModel requires sklearn to be installed to evaluate the model."
POTENTIAL_LOSS_COEF_MUST_BE_NON_NEGATIVE_ERROR = \
    "potential_loss_coef must be a non-negative float when include_potential is True, got {}."
ENERGY_LOSS_COEF_MUST_BE_NON_NEGATIVE_ERROR = \
    "energy_loss_coef must be a non-negative float, got {}."
LOADER_INCLUDE_POTENTIAL_MUST_BE_TRUE_ERROR = \
    "When include_potential is True, {} must have include_potential=True; got include_potential=False."
GRADIENT_CLIP_NORM_MUST_BE_POSITIVE_ERROR = \
    "gradient_clip_norm must be a positive float when gradient_clip is True, got {}."
SCALER_MUST_HAVE_MEAN_OR_CENTER_ERROR = \
    "scaler must have 'mean_' or 'center_' attribute; got type {}."
SCALER_MUST_HAVE_SCALE_ERROR = \
    "scaler must have 'scale_' attribute; got type {}."
INVALID_FEATURE_NAME_FOR_RESPONSE_KERNEL_COMPUTATION_ERROR = (
    "Invalid feature name for response kernel computation: {}. "
    "Valid feature names are: rho, grad_rho, grad_rho_norm, grad_rho_abs, "
    "grad_rho_mag, lap_rho, grad_rho_reduced, lap_rho_reduced."
)



# Warnings
PATIENCE_NOT_NONE_WHEN_VAL_LOADER_IS_NONE_WARNING = \
    "WARNING: patience should be None when val_loader is None, get {} instead."
INCLUDE_POTENTIAL_MUST_BE_NONE_FOR_POTENTIAL_WARNING = \
    "include_potential must be None when training a potential model (model_kind='potential'); ignored."
POTENTIAL_LOSS_COEF_MUST_BE_NONE_FOR_POTENTIAL_WARNING = \
    "potential_loss_coef must be None when training a potential model (model_kind='potential'); ignored."
ENERGY_LOSS_COEF_MUST_BE_NONE_FOR_POTENTIAL_WARNING = \
    "energy_loss_coef must be None when training a potential model (model_kind='potential'); ignored."
ENERGY_TARGET_MUST_BE_NONE_FOR_POTENTIAL_WARNING = \
    "energy_target must be None when training a potential model (model_kind='potential'); ignored."
ENERGY_TARGET_NOT_VALID_ERROR = \
    "energy_target must be 'total_energy' or 'energy_density', got {}."
TOTAL_ENERGY_TRUE_REQUIRED_ERROR = \
    "total_energy_true is required when energy_target='total_energy'. Use ExcDataManager.prepare_energy_dataloader() to create the loader."

# exclude_rho_from_nn warnings
EXCLUDE_RHO_REQUIRES_ENHANCEMENT_FACTOR_WARNING = \
    "exclude_rho_from_nn=True is intended for enhancement factor fitting; target_component is '{}', not enhancement_factor_x."
EXCLUDE_RHO_REQUIRES_REDUCED_FEATURES_WARNING = \
    "exclude_rho_from_nn=True expects reduced features (grad_rho_reduced, lap_rho_reduced); features_list has non-reduced: {}."



# Train() parameter types (extend Literal when adding optimizers/criteria)
OptimizerType    = Literal["Adam", "SGD"]
CriterionType    = Literal["L1Loss", "MSELoss", "L1norm", "L2norm"]
EnergyTargetType = Literal["total_energy", "energy_density"]

# Default batch sizes: potential = samples per batch; energy = configurations per batch
DEFAULT_BATCH_SIZE_POTENTIAL = 256
DEFAULT_BATCH_SIZE_ENERGY = 3

# Default loss coefficients for energy training (depend on energy_target)
DEFAULT_POTENTIAL_LOSS_COEF_TOTAL_ENERGY   = 1e2   # total_energy: strong potential gradient (dE/dρ→v_xc) matching
DEFAULT_POTENTIAL_LOSS_COEF_ENERGY_DENSITY = 1.0   # energy_density: potential helps gradient matching
DEFAULT_ENERGY_LOSS_COEF_TOTAL_ENERGY     = 1.0   # total_energy: integrated E is main target
DEFAULT_ENERGY_LOSS_COEF_ENERGY_DENSITY   = 1.0   # energy_density: grid-wise energy density

# Gradient clipping: set DEFAULT_GRADIENT_CLIP = True to enable by default; DEFAULT_GRADIENT_CLIP_NORM is used when enabled
DEFAULT_GRADIENT_CLIP = False
DEFAULT_GRADIENT_CLIP_NORM = 1.0


class TorchXCModel(XCBaseModel):
    """
    Unified Torch model for both potential and energy.
    """

    model : nn.Module = field(init=False)

    def __init__(
        self,
        *,
        model_kind        : ModelKind,
        features_list     : List[str],
        model_cls         : Optional[Type["nn.Module"]] = None,
        model_init_kwargs : Optional[Dict[str, Any]]    = None,
        model_instance    : Optional["nn.Module"]       = None,
        model_name        : Optional[str]               = None,
        model_dir         : str                         = "./models",
        model_config      : Optional[Dict[str, Any]]    = None,
        config_ext        : ConfigExt                   = "json",
        device            : Optional[str]               = "cpu",
        exclude_rho_from_nn: bool                       = False,
    ) -> None:

        if not TORCH_AVAILABLE:
            raise ImportError(TORCH_NOT_AVAILABLE_FOR_TORCHXCMODEL_ERROR)

        # When exclude_rho_from_nn=True: validate and override input_dim
        model_init_kwargs_to_use = model_init_kwargs
        if exclude_rho_from_nn and "rho" in features_list:
            non_rho = [f for f in features_list if f != "rho"]
            reduced_only = {"grad_rho_reduced", "lap_rho_reduced"}
            if not all(f in reduced_only for f in non_rho):
                warnings.warn(
                    EXCLUDE_RHO_REQUIRES_REDUCED_FEATURES_WARNING.format(
                        [f for f in non_rho if f not in reduced_only]
                    )
                )
            model_init_kwargs_to_use = dict(model_init_kwargs or {})
            model_init_kwargs_to_use["input_dim"] = len(features_list) - 1

        super().__init__(
            model_kind        = model_kind,
            features_list     = features_list,
            weights_ext       = "pth",
            config_ext        = config_ext,
            model_cls         = model_cls,
            model_init_kwargs = model_init_kwargs_to_use,
            model_instance    = model_instance,
            model_name        = model_name,
            model_dir         = model_dir,
            device            = device,
        )

        self.exclude_rho_from_nn = exclude_rho_from_nn

        if not isinstance(self.model, nn.Module):
            raise TypeError(MODEL_NOT_NNMODULE_ERROR.format(type(self.model)))

        if device is not None:
            self.model.to(device)
        self.model.double()  # use float64 for higher precision


    def train(
        self,
        train_loader        : DataLoaderType,
        val_loader          : Optional[DataLoaderType] = None,
        epochs              : int                      = 100,
        lr                  : float                    = 1e-3,
        batch_size          : Optional[int]            = None,
        optimizer           : OptimizerType            = "Adam",
        criterion           : CriterionType            = "L1Loss",
        patience            : Optional[int]            = 100,
        optimizer_kwargs    : Dict[str, Any]           = {},
        seed                : Optional[int]            = 42,
        include_potential   : Optional[bool]           = None,
        potential_loss_coef : Optional[float]          = None,
        energy_loss_coef    : Optional[float]          = None,
        gradient_clip       : Optional[bool]           = None,
        gradient_clip_norm  : Optional[float]          = None,
        energy_target       : EnergyTargetType         = "total_energy",
    ):
        """
        Train the model on the given data loader.

        Parameters
        ----------
        train_loader : VxcDataLoader | ExcDataLoader
            Training data; type must match self.model_kind (potential → Vxc, energy → Exc).
        val_loader : VxcDataLoader | ExcDataLoader, optional
            Validation data; type must match train_loader. If invalid, set to None with a warning.
        epochs : int
            Number of training epochs.
        lr : float
            Learning rate for the optimizer.
        batch_size : int, optional
            If None: potential → 256 (samples per batch), energy → 3 (configurations per batch).
        optimizer : OptimizerType
            "Adam" or "SGD".
        criterion : CriterionType
            Per-grid loss: ``L1Loss`` or ``L1norm`` → MAE; ``MSELoss`` or ``L2norm`` → MSE (default ``L1Loss``).
            Pair L1 variants with ``loss_norm_for_potential='L1norm'``; L2/MSE variants with ``'L2norm'``.
        patience : int, optional
            Early stopping patience (number of epochs without improvement). If val_loader is None, set to epochs+1.
        optimizer_kwargs : dict
            Extra keyword arguments for the optimizer constructor.
        seed : int, optional
            Random seed for shuffling and dropout.
        include_potential : bool, optional
            (Energy only.) Whether to include the potential-derived loss term (dE/dρ → v_xc/v_x/v_c).
            If None and model_kind == "energy", defaults to True. Must be None when model_kind == "potential".
        potential_loss_coef : float, optional
            (Energy only.) Coefficient for the potential-derived loss term in the total loss
            (total = energy_loss_coef * loss_energy + potential_loss_coef * loss_potential). Distinct from per-grid sample weights.
            If None: total_energy → 1e2, energy_density → 1.0. Must be None when model_kind == "potential".
            Must be >= 0 when include_potential is True.
        energy_loss_coef : float, optional
            (Energy only.) Coefficient for the energy density loss term in the total loss.
            If None: total_energy → 1.0, energy_density → 1.0. Must be None when model_kind == "potential".
            Must be >= 0.
        gradient_clip : bool, optional
            (Energy only.) Whether to apply gradient clipping. If None, uses DEFAULT_GRADIENT_CLIP.
        gradient_clip_norm : float, optional
            (Energy only.) Max norm for gradient clipping when gradient_clip is True.
            If None when gradient_clip is True, uses DEFAULT_GRADIENT_CLIP_NORM.
        energy_target : {"total_energy", "energy_density"}, optional
            (Energy only.) Target for the energy loss term.
            - "total_energy": minimize loss on total energy (integrated over grid).
            - "energy_density": minimize loss on grid-wise energy density.
            Default is "total_energy".

        Returns
        -------
        model, train_loss_list, val_loss_list
        """
        # When exclude_rho_from_nn: warn if target is not enhancement factor
        if self.exclude_rho_from_nn:
            target = getattr(train_loader, "target_component", None)
            if target != "enhancement_factor_x":
                warnings.warn(
                    EXCLUDE_RHO_REQUIRES_ENHANCEMENT_FACTOR_WARNING.format(target or "unknown")
                )

        # Type checks for energy-only args based on self.model_kind (not train_loader)
        if self.model_kind == "potential":
            if include_potential is not None:
                warnings.warn(INCLUDE_POTENTIAL_MUST_BE_NONE_FOR_POTENTIAL_WARNING)
            if potential_loss_coef is not None:
                warnings.warn(POTENTIAL_LOSS_COEF_MUST_BE_NONE_FOR_POTENTIAL_WARNING)
            if energy_loss_coef is not None:
                warnings.warn(ENERGY_LOSS_COEF_MUST_BE_NONE_FOR_POTENTIAL_WARNING)
            if energy_target != "total_energy":
                warnings.warn(ENERGY_TARGET_MUST_BE_NONE_FOR_POTENTIAL_WARNING)
        else:  # energy
            if include_potential is None:
                include_potential = True
            if potential_loss_coef is None:
                potential_loss_coef = (DEFAULT_POTENTIAL_LOSS_COEF_TOTAL_ENERGY if energy_target == "total_energy"
                                      else DEFAULT_POTENTIAL_LOSS_COEF_ENERGY_DENSITY)
            if energy_loss_coef is None:
                energy_loss_coef = (DEFAULT_ENERGY_LOSS_COEF_TOTAL_ENERGY if energy_target == "total_energy"
                                   else DEFAULT_ENERGY_LOSS_COEF_ENERGY_DENSITY)
            if include_potential and (not isinstance(potential_loss_coef, (int, float)) or potential_loss_coef < 0):
                raise ValueError(POTENTIAL_LOSS_COEF_MUST_BE_NON_NEGATIVE_ERROR.format(potential_loss_coef))
            if not isinstance(energy_loss_coef, (int, float)) or energy_loss_coef < 0:
                raise ValueError(ENERGY_LOSS_COEF_MUST_BE_NON_NEGATIVE_ERROR.format(energy_loss_coef))
            if energy_target not in ("total_energy", "energy_density"):
                raise ValueError(ENERGY_TARGET_NOT_VALID_ERROR.format(energy_target))
            # When include_potential is True, loaders must also have include_potential=True
            if include_potential:
                if isinstance(train_loader, ExcDataLoader) and not train_loader.include_potential:
                    raise ValueError(LOADER_INCLUDE_POTENTIAL_MUST_BE_TRUE_ERROR.format("train_loader"))
                if val_loader is not None and isinstance(val_loader, ExcDataLoader) and not val_loader.include_potential:
                    raise ValueError(LOADER_INCLUDE_POTENTIAL_MUST_BE_TRUE_ERROR.format("val_loader"))

        # Resolve batch_size default by self.model_kind
        if batch_size is None:
            batch_size = DEFAULT_BATCH_SIZE_POTENTIAL if self.model_kind == "potential" else DEFAULT_BATCH_SIZE_ENERGY

        # Resolve gradient_clip for energy training
        if self.model_kind == "energy":
            if gradient_clip is None:
                gradient_clip = DEFAULT_GRADIENT_CLIP
            if gradient_clip and gradient_clip_norm is None:
                gradient_clip_norm = DEFAULT_GRADIENT_CLIP_NORM
            if gradient_clip and (not isinstance(gradient_clip_norm, (int, float)) or gradient_clip_norm <= 0):
                raise ValueError(GRADIENT_CLIP_NORM_MUST_BE_POSITIVE_ERROR.format(gradient_clip_norm))

        # check the validation data loader (must match model_kind: potential -> Vx
        # c, energy -> Exc)
        if val_loader is not None:
            try:
                if self.model_kind == "potential":
                    self._check_vxc_data_loader(val_loader)
                elif self.model_kind == "energy":
                    self._check_exc_data_loader(val_loader)
                else:
                    raise ValueError(DATA_LOADER_TYPE_NOT_VALID_ERROR.format(self.model_kind))
            except ValueError as e:
                print(f"WARNING: {e}")
                print("\tThe validation loader will set to be None for this run.")
                val_loader = None

        # set patience to be larger than epochs if validation loader is not provided
        if val_loader is None and patience is not None:
            print(PATIENCE_NOT_NONE_WHEN_VAL_LOADER_IS_NONE_WARNING.format(patience))
            print("\tThe patience will be set to be larger than the number of epochs for this run.")
            patience = epochs + 1

        # check parameters
        self.check_epochs(epochs)
        self.check_lr(lr)
        self.check_patience(patience)

        # check and get the optimizer (skip training if model has no learnable parameters)
        params = list(self.model.parameters())
        if len(params) == 0:
            warnings.warn(
                "Model has no learnable parameters. Skipping training. "
                "If using FixedPBEEnhancementFactor, ensure nn_models is reloaded "
                "(importlib.reload(nn_models)) or restart the kernel after adding the scale parameter.",
                UserWarning,
                stacklevel=2,
            )

        if optimizer == "Adam":
            optimizer = torch.optim.Adam(params, lr=lr, **optimizer_kwargs)
        elif optimizer == "SGD":
            optimizer = torch.optim.SGD(params, lr=lr, **optimizer_kwargs)
        else:
            raise ValueError(OPTIMIZER_NOT_VALID_ERROR.format(optimizer))
        
        # check and get the criterion (L1norm/L2norm mirror naming in data_manager weights)
        if criterion in ("L1Loss", "L1norm"):
            criterion = torch.nn.L1Loss(reduction='none')
        elif criterion in ("MSELoss", "L2norm"):
            criterion = torch.nn.MSELoss(reduction='none')
        else:
            raise ValueError(CRITERION_NOT_VALID_ERROR.format(criterion))
        
        # train the model
        if self.model_kind == "potential":
            return self._train_potential(
                train_loader = train_loader,
                val_loader   = val_loader,
                epochs       = epochs,
                batch_size   = batch_size,
                patience     = patience,
                optimizer    = optimizer,
                criterion    = criterion,
                seed         = seed,
            )
        elif self.model_kind == "energy":
            return self._train_energy(
                train_loader = train_loader,
                val_loader   = val_loader,
                epochs       = epochs,
                batch_size   = batch_size,
                patience     = patience,
                optimizer    = optimizer,
                criterion    = criterion,
                seed         = seed,
                # energy-only parameters
                include_potential   = include_potential,
                potential_loss_coef = potential_loss_coef,
                energy_loss_coef    = energy_loss_coef,
                gradient_clip       = gradient_clip,
                gradient_clip_norm  = gradient_clip_norm,
                energy_target       = energy_target,
            )
        else:
            raise ValueError(f"Invalid model kind: {self.model_kind}")


    def _train_potential(
        self,
        train_loader : VxcDataLoader,
        val_loader   : Optional[VxcDataLoader],
        epochs       : int,
        batch_size   : int,
        patience     : Optional[int],
        optimizer    : torch.optim.Optimizer,
        criterion    : torch.nn.Module,
        seed         : Optional[int] = None,
    ):
        # check and transform the training data loader
        self._check_vxc_data_loader(train_loader)
        train_loader = self.transform_to_torch_loader(train_loader, batch_size, shuffle=True , seed=seed)
        val_loader   = self.transform_to_torch_loader(val_loader  , batch_size, shuffle=False, seed=None) if val_loader is not None else None

        # train the model
        best_val_loss = float('inf')
        patience_counter = 0
        best_model_state = None

        # loss list
        train_loss_list = []
        val_loss_list = []

        for epoch in range(epochs):
            # training
            self.model.train()
            train_loss = 0.0

            batch_count = 0
            for batch_data in train_loader:
                train_X, train_y, train_weights = batch_data
                train_weights = train_weights.squeeze()

                optimizer.zero_grad()
                nn_indices = self._get_nn_feature_indices()
                train_X_for_nn = train_X[:, nn_indices]
                outputs = self.model(train_X_for_nn)
                per_sample_loss = criterion(outputs, train_y).squeeze()

                assert train_weights.shape == per_sample_loss.shape, \
                    WEIGHTS_NOT_SAME_AS_SAMPLES_ERROR.format(train_weights.shape, per_sample_loss.shape)

                loss = (per_sample_loss * train_weights).mean()
                
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
                batch_count += 1
            train_loss = train_loss / max(batch_count, 1)
            train_loss_list.append(train_loss)
            

            # validation (in batches to avoid loading full val set into memory)
            self.model.eval()
            val_loss = 0.0
            if val_loader is not None:
                device = self.device if self.device else "cpu"
                val_loss_sum = 0.0
                val_n = 0
                with torch.no_grad():
                    for batch_data in val_loader:
                        X_batch, y_batch, w_batch = batch_data
                        X_batch = X_batch.to(device)
                        y_batch = y_batch.to(device)
                        w_batch = w_batch.squeeze().to(device)
                        nn_indices = self._get_nn_feature_indices()
                        X_batch_for_nn = X_batch[:, nn_indices]
                        y_pred = self.model(X_batch_for_nn)
                        per_sample_loss = criterion(y_pred, y_batch).squeeze()
                        val_loss_sum += (per_sample_loss * w_batch).sum().item()
                        val_n += per_sample_loss.numel()
                val_loss = val_loss_sum / max(val_n, 1)
                val_loss_list.append(val_loss)
            else:
                val_loss_list.append(-1.0)

            # early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
                best_model_state = self.model.state_dict().copy()
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch + 1}")
                    break

            if (epoch + 1) % 1 == 0:
                print(f"Epoch {epoch+1:>4d}/{epochs}, Train Loss: {train_loss:.6e}, Val Loss: {val_loss:.6e}")

        # load best model
        if best_model_state is not None:
            self.model.load_state_dict(best_model_state)

        return self.model, train_loss_list, val_loss_list


    def _train_energy(
        self,
        train_loader        : ExcDataLoader,
        val_loader          : Optional[ExcDataLoader],
        epochs              : int,
        batch_size          : int,
        patience            : Optional[int],
        optimizer           : Optional[torch.optim.Optimizer],
        criterion           : Optional[torch.nn.Module],
        seed                : Optional[int]    = None,
        include_potential   : bool             = True,
        potential_loss_coef : float            = DEFAULT_POTENTIAL_LOSS_COEF_TOTAL_ENERGY,
        energy_loss_coef    : float            = DEFAULT_ENERGY_LOSS_COEF_TOTAL_ENERGY,
        gradient_clip       : bool             = False,
        gradient_clip_norm  : Optional[float]  = None,
        energy_target       : EnergyTargetType = "total_energy",
    ):        
        # check and transform the training data loader
        self._check_exc_data_loader(train_loader)
        device = self.device if self.device else "cpu"
        shared_derivative_matrix = train_loader.shared_derivative_matrix
        computation_config = train_loader.get_computation_config()
        train_loader = self.transform_to_torch_loader(train_loader, batch_size, shuffle=True, seed=seed)

        # train the model
        best_val_loss    = float('inf')
        patience_counter = 0
        best_model_state = None

        # loss list
        train_loss_list = []
        val_loss_list = []

        for epoch in range(epochs):
            # training
            self.model.train()
            train_loss = 0.0
            train_energy_loss = 0.0
            train_potential_loss = 0.0
            config_count = 0  # Count total configurations, not batches
            for batch_data in train_loader:
                
                # zero the gradients
                optimizer.zero_grad()
                loss = 0.0 
                
                # compute the loss for each configuration
                for _, config_info_dict in enumerate(batch_data):
                    # Get the configuration information and apply grid cutoff
                    max_element_cutoff_index = config_info_dict["max_element_cutoff_index"]
                    quadrature_point_number  = config_info_dict["quadrature_point_number"]
                    max_grid_cutoff_index    = max_element_cutoff_index * quadrature_point_number

                    # Get the configuration data
                    quadrature_nodes   = config_info_dict["quadrature_nodes"][:max_grid_cutoff_index]
                    quadrature_weights = config_info_dict["quadrature_weights"][:max_grid_cutoff_index]
                    features_physical  = config_info_dict["features_physical"][:max_grid_cutoff_index]
                    y_energy_true      = config_info_dict["energy_density_transformed"][:max_grid_cutoff_index]
                    y_potential_true   = config_info_dict["potential_transformed"][:max_grid_cutoff_index]
                    energy_weights     = config_info_dict["weights_for_energy_density"][:max_grid_cutoff_index]
                    potential_weights  = config_info_dict["weights_for_potential"][:max_grid_cutoff_index]
                    derivative_matrix  = config_info_dict["derivative_matrix"]
                    if derivative_matrix is None:
                        derivative_matrix = shared_derivative_matrix
                    derivative_matrix = derivative_matrix[:max_element_cutoff_index]

                    # Prepare total_energy_true for total_energy loss (convert to tensor and move to device)
                    E_true_raw = config_info_dict.get("total_energy_true")
                    if energy_target == "total_energy":
                        if E_true_raw is None:
                            raise ValueError(TOTAL_ENERGY_TRUE_REQUIRED_ERROR)
                        E_true = torch.as_tensor(E_true_raw, dtype=torch.float64, device=device)
                    else:
                        E_true = None

                    # Move tensors to model device (collate_fn creates on CPU)                    
                    features_physical  = features_physical.to(device)
                    y_energy_true      = y_energy_true.to(device)
                    energy_weights     = energy_weights.to(device)
                    quadrature_nodes   = quadrature_nodes.to(device)
                    quadrature_weights = quadrature_weights.to(device)
                    derivative_matrix  = derivative_matrix.to(device)
                    if include_potential:
                        # move the potential to the model device
                        y_potential_true  = y_potential_true.to(device)
                        potential_weights = potential_weights.to(device)

                        # Ensure features_physical requires gradients for potential computation
                        # This must be done before calling model() so the computation graph includes features
                        features_physical = features_physical.requires_grad_(True)

                    y_energy_physical, y_energy_transformed = self.compute_energy_density_from_features(
                        model             = self.model,
                        features_physical = features_physical,
                        config            = computation_config,
                        debug_context     = "Training",
                        exclude_rho_from_nn = self.exclude_rho_from_nn,
                    )

                    # Compute the energy loss for the current configuration
                    assert energy_weights.shape == y_energy_transformed.shape, \
                        WEIGHTS_NOT_SAME_AS_SAMPLES_ERROR.format(energy_weights.shape, y_energy_transformed.shape)
                    if energy_target == "total_energy":
                        # Loss on integrated total energy; E_true precomputed in prepare_energy_dataloader
                        r_sq_w = 4 * np.pi * (quadrature_nodes ** 2) * quadrature_weights
                        if y_energy_physical.dim() > 1:
                            r_sq_w = r_sq_w.unsqueeze(1)
                        E_pred = (r_sq_w * y_energy_physical).sum(dim=0).reshape(1, -1)
                        E_true = E_true.reshape(1, -1)
                        energy_loss = criterion(E_pred, E_true).mean()
                    else:
                        # Loss on grid-wise energy density (energy_target == "energy_density")
                        energy_per_sample_loss = criterion(y_energy_transformed, y_energy_true).squeeze()
                        energy_loss = (energy_per_sample_loss * energy_weights).mean()

                    # Then, compute the predicted potential (accepts physical features)
                    if include_potential:
                        y_potential_pred = self.compute_potential_from_energy_density(
                            y_energy_physical   = y_energy_physical,
                            features_physical   = features_physical,
                            derivative_matrix   = derivative_matrix,
                            quadrature_nodes    = quadrature_nodes,
                            config              = computation_config,
                            requires_grad       = True,  # Training: need gradients for backpropagation
                            transform_potential = True,  # Training: transform to match targets
                        )
                        
                        # Debug: Check y_potential_pred and y_potential_true for NaN/Inf
                        self._debug_check_nan_inf(y_potential_pred, "y_potential_pred", context="Training", include_stats=True)
                        self._debug_check_nan_inf(y_potential_true, "y_potential_true", context="Training")

                        # Compute the potential loss for the current configuration, if needed
                        assert potential_weights.shape == y_potential_pred.shape, \
                            WEIGHTS_NOT_SAME_AS_SAMPLES_ERROR.format(potential_weights.shape, y_potential_pred.shape)
                        potential_per_sample_loss = criterion(y_potential_pred, y_potential_true).squeeze()
                        potential_loss = (potential_per_sample_loss * potential_weights).mean()
                        
                        # Debug: Check potential_loss for NaN/Inf
                        self._debug_check_nan_inf(
                            potential_loss, "potential_loss",
                            context="Training",
                            extra=f"y_potential_pred shape={y_potential_pred.shape}, y_potential_true shape={y_potential_true.shape}, "
                                  f"y_potential_pred range=[{y_potential_pred.min()}, {y_potential_pred.max()}], "
                                  f"y_potential_true range=[{y_potential_true.min()}, {y_potential_true.max()}]",
                        )
                    else:
                        potential_loss = 0.0
                    
                    # update the total train loss
                    batch_loss = energy_loss_coef * energy_loss + potential_loss_coef * potential_loss
                    
                    # Debug: Check batch_loss for NaN/Inf
                    self._debug_check_nan_inf(
                        batch_loss, "batch_loss",
                        context="Training",
                        extra=f"energy_loss={energy_loss}, potential_loss={potential_loss}, energy_loss_coef={energy_loss_coef}, potential_loss_coef={potential_loss_coef}",
                    )
                    
                    loss += batch_loss
                    # Accumulate energy and potential losses separately
                    train_energy_loss += energy_loss.item()
                    if include_potential:
                        train_potential_loss += potential_loss.item()
                    config_count += 1  # Count configurations, not batches

                # Backward pass after accumulating all losses in the batch
                loss.backward()
                
                if gradient_clip and gradient_clip_norm is not None:
                    torch.nn.utils.clip_grad_norm_(self.model.parameters(), gradient_clip_norm)
                
                # Debug: Check gradients for NaN/Inf before optimizer step
                has_nan_grad = False
                has_inf_grad = False
                for name, param in self.model.named_parameters():
                    if param.grad is not None and self._debug_check_nan_inf(param.grad, f"Gradient for '{name}'", context="Training"):
                        has_nan_grad = True
                        has_inf_grad = True
                if has_nan_grad or has_inf_grad:
                    print(f"DEBUG Training: Skipping optimizer step due to NaN/Inf gradients")
                    # Zero gradients to prevent accumulation
                    optimizer.zero_grad()
                    continue
                
                optimizer.step()
                
                # Debug: Check total loss before item()
                self._debug_check_nan_inf(loss, "total loss before item()", context="Training", extra=f"loss={loss}")
                
                # Accumulate loss per configuration (not per batch)
                # loss is the sum of all configurations in this batch
                train_loss += loss.item()

            # Average over total number of configurations
            train_loss = train_loss / max(config_count, 1)
            train_energy_loss = train_energy_loss / max(config_count, 1)
            if include_potential:
                train_potential_loss = train_potential_loss / max(config_count, 1)
            train_loss_list.append(train_loss)

            # Validation
            self.model.eval()
            val_loss           = 0.0
            val_energy_loss    = 0.0
            val_potential_loss = 0.0
            batch_count        = 0
            if val_loader is not None:
                for _, exc_configuration in enumerate(val_loader):
                    # Get the configuration information and apply grid cutoff
                    max_element_cutoff_index = exc_configuration.max_element_cutoff_index
                    quadrature_point_number  = exc_configuration.quadrature_point_number
                    max_grid_cutoff_index    = max_element_cutoff_index * quadrature_point_number

                    # Get the configuration data
                    quadrature_nodes   = torch.DoubleTensor(exc_configuration.quadrature_nodes[:max_grid_cutoff_index]).to(device)
                    quadrature_weights = torch.DoubleTensor(exc_configuration.quadrature_weights[:max_grid_cutoff_index]).to(device)
                    features_physical  = torch.DoubleTensor(exc_configuration.features_physical[:max_grid_cutoff_index]).to(device)
                    y_energy_true      = torch.DoubleTensor(exc_configuration.energy_density_transformed[:max_grid_cutoff_index]).to(device)
                    energy_weights     = torch.DoubleTensor(exc_configuration.weights_for_energy_density[:max_grid_cutoff_index]).to(device)
                    derivative_matrix  = exc_configuration.derivative_matrix
                    if derivative_matrix is None:
                        derivative_matrix = shared_derivative_matrix
                    derivative_matrix = torch.DoubleTensor(derivative_matrix[:max_element_cutoff_index]).to(device)

                    # Prepare total_energy_true for total_energy loss (convert to tensor and move to device)
                    E_true_raw = exc_configuration.total_energy_true
                    if energy_target == "total_energy":
                        if E_true_raw is None:
                            raise ValueError(TOTAL_ENERGY_TRUE_REQUIRED_ERROR)
                        E_true = torch.as_tensor(E_true_raw, dtype=torch.float64, device=device)
                    else:
                        E_true = None

                    if include_potential:
                        y_potential_true  = torch.DoubleTensor(exc_configuration.potential_transformed[:max_grid_cutoff_index]).to(device)
                        potential_weights = torch.DoubleTensor(exc_configuration.weights_for_potential[:max_grid_cutoff_index]).to(device)
                        features_physical = features_physical.requires_grad_(True)

                    y_energy_physical, y_energy_transformed = self.compute_energy_density_from_features(
                        model             = self.model,
                        features_physical = features_physical,
                        config            = computation_config,
                        debug_context     = "Validation",
                        exclude_rho_from_nn = self.exclude_rho_from_nn,
                    )

                    # Compute the energy loss for the current configuration
                    assert energy_weights.shape == y_energy_transformed.shape, \
                        WEIGHTS_NOT_SAME_AS_SAMPLES_ERROR.format(energy_weights.shape, y_energy_transformed.shape)
                    if energy_target == "total_energy":
                        r_sq_w = 4 * np.pi * (quadrature_nodes ** 2) * quadrature_weights
                        if y_energy_physical.dim() > 1:
                            r_sq_w = r_sq_w.unsqueeze(1)
                        E_pred = (r_sq_w * y_energy_physical).sum(dim=0).reshape(1, -1)
                        E_true = E_true.reshape(1, -1)
                        energy_loss = criterion(E_pred, E_true).mean()
                    else:
                        energy_per_sample_loss = criterion(y_energy_transformed, y_energy_true).squeeze()
                        energy_loss = (energy_per_sample_loss * energy_weights).mean()

                    # Then, compute the predicted potential, which is more complicated
                    if include_potential:
                        y_potential_pred = self.compute_potential_from_energy_density(
                            y_energy_physical   = y_energy_physical,
                            features_physical   = features_physical,
                            derivative_matrix   = derivative_matrix,
                            quadrature_nodes    = quadrature_nodes,
                            config              = computation_config,
                            requires_grad       = False,  # Validation: no gradients needed for backprop, but needed for computation
                            transform_potential = True,   # Validation: must match y_potential_true which is in transformed space
                        )
                        
                        # Debug: Check y_potential_pred and y_potential_true for NaN/Inf
                        self._debug_check_nan_inf(y_potential_pred, "y_potential_pred", context="Validation", include_stats=True)
                        self._debug_check_nan_inf(y_potential_true, "y_potential_true", context="Validation")

                        # Compute the potential loss for the current configuration
                        assert potential_weights.shape == y_potential_pred.shape, \
                            WEIGHTS_NOT_SAME_AS_SAMPLES_ERROR.format(potential_weights.shape, y_potential_pred.shape)
                        potential_per_sample_loss = criterion(y_potential_pred, y_potential_true).squeeze()
                        potential_loss = (potential_per_sample_loss * potential_weights).mean()
                        
                        # Debug: Check potential_loss for NaN/Inf
                        self._debug_check_nan_inf(
                            potential_loss, "potential_loss",
                            context="Validation",
                            extra=f"y_potential_pred shape={y_potential_pred.shape}, y_potential_true shape={y_potential_true.shape}, "
                                  f"y_potential_pred range=[{y_potential_pred.min()}, {y_potential_pred.max()}], "
                                  f"y_potential_true range=[{y_potential_true.min()}, {y_potential_true.max()}]",
                        )
                    else:
                        potential_loss = 0.0
                    
                    # Use no_grad only for loss accumulation (not for computation)
                    with torch.no_grad():
                        val_loss_batch = energy_loss_coef * energy_loss + potential_loss_coef * potential_loss
                        
                        # Debug: Check validation loss batch for NaN/Inf
                        self._debug_check_nan_inf(
                            val_loss_batch, "val_loss_batch",
                            context="Validation",
                            extra=f"energy_loss={energy_loss}, potential_loss={potential_loss}, energy_loss_coef={energy_loss_coef}, potential_loss_coef={potential_loss_coef}",
                        )
                        
                        val_loss += val_loss_batch.item()
                        # Accumulate energy and potential losses separately
                        val_energy_loss += energy_loss.item()
                        if include_potential:
                            val_potential_loss += potential_loss.item()
                    batch_count += 1

            if batch_count > 0:
                val_loss = val_loss / batch_count
                val_energy_loss = val_energy_loss / batch_count
                if include_potential:
                    val_potential_loss = val_potential_loss / batch_count
            val_loss_list.append(val_loss)  # val_loss is already a float, no need for .item()

            # early stopping
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
                best_model_state = self.model.state_dict().copy()
            else:
                patience_counter += 1
                if patience is not None and patience_counter >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break

            if (epoch + 1) % 1 == 0:
                if include_potential:
                    print(f"Epoch {epoch+1:>4d}/{epochs}, Train Loss: {train_loss:.6e} (Energy: {train_energy_loss:.6e}, Potential: {train_potential_loss:.6e}), Val Loss: {val_loss:.6e} (Energy: {val_energy_loss:.6e}, Potential: {val_potential_loss:.6e})")
                else:
                    print(f"Epoch {epoch+1:>4d}/{epochs}, Train Loss: {train_loss:.6e}, Val Loss: {val_loss:.6e}")

        # load best model
        if best_model_state is not None:
            self.model.load_state_dict(best_model_state)

        return self.model, train_loss_list, val_loss_list



    @classmethod
    def compute_energy_density_from_features(
        cls,
        model             : torch.nn.Module,
        features_physical : torch.Tensor,
        config            : ComputationConfig,
        debug_context     : Optional[str] = None,
        exclude_rho_from_nn: bool = False,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute energy density from features.

        Applies differentiable transform to features, feeds to model, then if target
        is enhancement factor, multiplies by LDA energy density.

        Parameters
        ----------
        model : torch.nn.Module
            Neural network model for energy density / enhancement factor prediction.
        features_physical : torch.Tensor
            Physical features, shape (n_grid, n_features).
        config : ComputationConfig
            Configuration with transform and target_component.
        debug_context : str, optional
            If provided (e.g. "Training", "Validation"), run NaN/Inf checks on inputs and outputs.

        Returns
        -------
        Tuple[torch.Tensor, torch.Tensor]
            (y_energy_physical, y_energy_transformed). Use y_energy_transformed for energy loss,
            y_energy_physical for compute_potential_from_energy_density (with features_physical).
        """
        features_list           = config.features_list
        scale_features          = config.scale_features
        scaler_for_features     = config.scaler_for_features
        use_symlog_for_features = config.use_symlog_for_features
        linthresh_for_features  = config.linthresh_for_features
        scale_energy            = config.scale_energy
        scaler_for_energy       = config.scaler_for_energy
        use_symlog_for_energy   = config.use_symlog_for_energy
        linthresh_for_energy    = config.linthresh_for_energy
        target_component        = config.target_component_for_energy

        # Debug: check inputs
        if debug_context:
            cls._debug_check_nan_inf(features_physical, "features_physical before model", context=debug_context, include_stats=True)
            model_has_nan = False
            for name, param in model.named_parameters():
                if cls._debug_check_nan_inf(param, f"Model parameter '{name}' before forward pass", context=debug_context):
                    model_has_nan = True
            if model_has_nan:
                print(f"DEBUG {debug_context}: Model has NaN/Inf parameters, this will cause NaN outputs!")

        # 1. Transform features (physical -> transformed)
        features_transformed = features_physical.clone()
        if use_symlog_for_features and linthresh_for_features is not None:
            features_transformed = cls._torch_symlog(features_transformed, linthresh_for_features)
        if scale_features and scaler_for_features is not None:
            features_transformed = cls._torch_apply_scaler(features_transformed, scaler_for_features)

        # 2. Model forward (exclude rho from NN input when exclude_rho_from_nn=True)
        if exclude_rho_from_nn and "rho" in features_list:
            nn_indices = [i for i, f in enumerate(features_list) if f != "rho"]
            features_for_nn = features_transformed[:, nn_indices]
        else:
            features_for_nn = features_transformed
        y_energy_pred = model(features_for_nn)

        # Debug: check model output
        if debug_context:
            if cls._debug_check_nan_inf(y_energy_pred, "y_energy_pred immediately after model", context=debug_context, include_stats=True):
                for name, param in model.named_parameters():
                    cls._debug_check_nan_inf(param, f"Model parameter '{name}' after forward pass", context=debug_context)

        # 3. Inverse transform energy density (transformed -> physical)
        y_energy_physical = y_energy_pred
        if scale_energy and scaler_for_energy is not None:
            y_energy_physical = cls._torch_inverse_scaler(y_energy_physical, scaler_for_energy)
        if use_symlog_for_energy and linthresh_for_energy is not None:
            y_energy_physical = cls._torch_symexp(y_energy_physical, linthresh_for_energy)


        # 4. If target is enhancement factor, multiply by LDA energy density
        if target_component == "enhancement_factor_x":
            rho = features_physical[:, features_list.index("rho")].reshape(-1, 1)
            lda_energy_density = cls._torch_lda_exchange_energy_density(rho, use_rho_43=True)
            y_energy_physical = y_energy_physical * lda_energy_density

        return y_energy_physical, y_energy_pred
        

    @classmethod
    def compute_potential_from_energy_density(
        cls,
        y_energy_physical   : torch.Tensor,
        features_physical   : torch.Tensor,
        derivative_matrix   : torch.Tensor,
        quadrature_nodes    : torch.Tensor,
        config              : ComputationConfig,
        requires_grad       : bool = True,
        transform_potential : bool = True,
    ) -> torch.Tensor:
        """
        Compute XC potential from energy density prediction using automatic differentiation
        or numerical differentiation.
        
        This method computes the functional derivative of the energy density with respect
        to the density (and its derivatives) to obtain the XC potential.
        
        For LDA: v_xc = δE_xc/δρ
        For GGA: v_xc = δE_xc/δρ - ∇·(δE_xc/δ(∇ρ))
        For meta-GGA: v_xc = δE_xc/δρ - ∇·(δE_xc/δ(∇ρ)) + ∇²(δE_xc/δ(∇²ρ)) + ...
        
        Parameters
        ----------
        y_energy_physical : torch.Tensor
            Predicted energy density values in physical space, shape (n_grid, n_targets) or (n_grid,).
            Should already include LDA multiplication if target is enhancement_factor_x.
        features_physical : torch.Tensor
            Input features in physical space, shape (n_grid, n_features).
            Contains density and its derivatives (rho, grad_rho, lap_rho, etc.).
        derivative_matrix : torch.Tensor
            Derivative matrix for computing gradients in spherical coordinates.
            Shape (n_elements, n_quad_per_element, n_quad_per_element) or (n_grid, n_grid).
            Used for computing ∇·(δE_xc/δ(∇ρ)) terms.
        quadrature_nodes : torch.Tensor
            Quadrature node positions (radial coordinates), shape (n_grid,).
            Used for spherical coordinate transformations.
        quadrature_weights : torch.Tensor
            Quadrature weights for numerical integration, shape (n_grid,).
            Used for computing functional derivatives.
        config : ComputationConfig
            Configuration object containing transformation and scaling parameters.
            Can be obtained from ExcDataLoader.get_computation_config() or
            MLXCCalculator.get_computation_config().
        requires_grad : bool, default=True
            Whether the returned potential tensor should require gradients.
            Set to True during training (to enable backpropagation through potential computation).
            Set to False during inference (to save memory and computation).
        transform_potential : bool, default=True
            Whether to transform the potential back to transformed space.
            Set to True during training (to match transformed target values).
            Set to False during inference (to return physical potential values).
        
        Returns
        -------
        torch.Tensor
            Predicted XC potential values, shape (n_grid, n_targets) or (n_grid,).
            Same shape as y_energy_pred.
            The returned tensor's requires_grad attribute is set according to requires_grad parameter.
        
        Notes
        -----
        This method implements the functional derivative computation for XC potentials.
        The computation involves:
        1. Extracting density and its derivatives from features
        2. Computing partial derivatives of energy density with respect to density variables
        3. Applying derivative operators for gradient terms: in spherical coordinates,
           -∇·(∂e/∂∇ρ) becomes -(1/r²)d/dr[r²·∂e/∂|∇ρ|]
        4. Transforming back from transformed space if needed
        
        The method supports both automatic differentiation (using torch.autograd) and
        numerical differentiation (finite differences).
        
        Examples
        --------
        >>> # Training: requires gradients for backpropagation
        >>> config = exc_data_loader.get_computation_config()
        >>> y_potential = TorchXCModel.compute_potential_from_energy_density(
        ...     y_energy_physical   = y_energy_physical,
        ...     features_physical  = features_physical,
        ...     derivative_matrix  = derivative_matrix,
        ...     quadrature_nodes   = quadrature_nodes,
        ...     config             = config,
        ...     requires_grad      = True,
        ... )
        
        >>> # Inference: no gradients needed
        >>> config = mlxc_calculator.get_computation_config()
        >>> y_potential = TorchXCModel.compute_potential_from_energy_density(
        ...     y_energy_physical   = y_energy_physical,
        ...     features_physical  = features_physical,
        ...     derivative_matrix  = derivative_matrix,
        ...     quadrature_nodes   = quadrature_nodes,
        ...     config             = config,
        ...     requires_grad      = False,
        ... )
        """

        # Extract configuration parameters
        features_list            = config.features_list
        target_component         = config.target_component_for_potential

        # optional parameters for scaling and symlog transformation
        scale_potential          = config.scale_potential
        scaler_for_potential     = config.scaler_for_potential
        use_symlog_for_potential = config.use_symlog_for_potential
        linthresh_for_potential  = config.linthresh_for_potential

        rho_physical = features_physical[:, features_list.index("rho")]

        # y_energy_physical is already in physical space (from compute_energy_density_from_features)
        y_energy_physical_pred = y_energy_physical

        # ===========================================================================================
        # Step 1: Compute pointwise gradients: ∂y_energy_physical_pred[i]/∂features_physical[i]
        #   Since we receive physical features and physical energy density, we compute gradients
        #   directly w.r.t. features_physical.
        # ===========================================================================================
        # Get dimensions
        n_grid = quadrature_nodes.shape[0]
        n_features = len(features_list)
        
        # Handle multiple targets: sum over targets for each grid point
        if y_energy_physical_pred.ndim > 1:
            energy_per_point = y_energy_physical_pred.sum(dim=1)  # (n_grid,)
        else:
            energy_per_point = y_energy_physical_pred  # (n_grid,)
        
        # Debug: Check for NaN/Inf in inputs
        cls._debug_check_nan_inf(y_energy_physical,      "y_energy_physical")
        cls._debug_check_nan_inf(y_energy_physical_pred, "y_energy_physical_pred")
        cls._debug_check_nan_inf(energy_per_point,       "energy_per_point")
        cls._debug_check_nan_inf(features_physical,      "features_physical")
        
        # Compute pointwise gradients: ∂energy_per_point[i]/∂features_physical[i]
        dE_dfeatures_physical = torch.autograd.grad(
            outputs      = energy_per_point,
            inputs       = features_physical,
            grad_outputs = torch.ones_like(energy_per_point),
            create_graph = requires_grad,
            retain_graph = requires_grad,
            only_inputs  = True,
        )[0]
        
        # Debug: Check for NaN/Inf in gradients
        if cls._debug_check_nan_inf(dE_dfeatures_physical, "dE_dfeatures_physical"):
            print(f"DEBUG: NaN locations: {torch.isnan(dE_dfeatures_physical).nonzero()}")
        
        # Check shape: (n_grid, n_features)
        assert dE_dfeatures_physical.shape == (n_grid, n_features), \
            f"Gradient shape mismatch: expected ({n_grid}, {n_features}), got {dE_dfeatures_physical.shape}"

        # ===========================================================================================
        # Step 2: Compute XC potential contributions from different feature channels
        #   Δv_xc = ∂Δe_xc/∂ρ - sign(∇ρ)·(1/r²)d/dr[r²·∂Δe_xc/∂|∇ρ|] + (1/r²)d/dr[r²·d/dr[∂Δe_xc/∂(∇²ρ)]]
        #   where (1/r²)d/dr[r²·f] is the radial divergence ∇·(f r̂) in spherical coordinates
        # ===========================================================================================
        # Initialize potential contributions
        vxc_contributions = torch.zeros(n_grid, dtype=features_physical.dtype, device=features_physical.device)

        # Process each feature channel
        for feat_idx, feat_name in enumerate(features_list):
            dE_dfeat = dE_dfeatures_physical[:, feat_idx]  # (n_grid,)
            
            # Debug: Check dE_dfeat for NaN/Inf
            cls._debug_check_nan_inf(dE_dfeat, f"dE_dfeat for {feat_name}")
            
            if feat_name == "rho":
                # Direct contribution: ∂Δe_xc/∂ρ
                vxc_contributions += dE_dfeat

                
            elif feat_name in ["grad_rho_norm", "grad_rho_abs", "grad_rho_mag"]:
                # Contribution: -sign(∇ρ)·(1/r²)d/dr[r²·∂Δe_xc/∂|∇ρ|]
                # sign(∇ρ) = -1 for |∇ρ|, so -sign(∇ρ)·(...) = +(...) → vxc += div_term
                div_term = cls._compute_radial_divergence(
                    dE_dfeat, derivative_matrix, quadrature_nodes
                )
                cls._debug_check_nan_inf(div_term, f"div_term for {feat_name}")
                vxc_contributions += div_term
                
            elif feat_name == "grad_rho":
                # Contribution: -d/dr[∂Δe_xc/∂∇ρ]
                div_term = cls._compute_radial_derivative(
                    dE_dfeat, derivative_matrix, quadrature_nodes
                )
                cls._debug_check_nan_inf(div_term, f"div_term for {feat_name}")
                vxc_contributions -= div_term
                
            elif feat_name == "lap_rho":
                # Contribution: (1/r²) * d/dr[r² * d/dr[∂Δe_xc/∂(∇²ρ)]]
                laplacian_contribution = cls._compute_radial_laplacian(
                    dE_dfeat, derivative_matrix, quadrature_nodes
                )
                # Debug: Check laplacian result
                cls._debug_check_nan_inf(laplacian_contribution, f"laplacian_contribution for {feat_name}")
                vxc_contributions += laplacian_contribution

            elif feat_name == "grad_rho_reduced":
                # Eq. 26 (appendix 3.1): s-contribution = -(4/3)(∂Δe/∂s)(s/ρ) - (1/r²)d/dr[r²·(∂Δe/∂s)(s/|∇ρ|)sign(∇ρ)]
                # s = |∇ρ|/(2·kf·ρ^(4/3)) => s/|∇ρ| = 1/(2·kf·ρ^(4/3))
                rho_idx = features_list.index("rho")
                rho = features_physical[:, rho_idx]
                s = features_physical[:, feat_idx]
                kf = (3 * np.pi**2) ** (1 / 3)
                rho_safe = torch.clamp(rho, min=1e-30)
                rho_43_safe = rho_safe ** (4 / 3)
                
                # Direct term: -(4/3)(∂Δe/∂s)(s/ρ)
                direct_term = - (4 / 3) * dE_dfeat * s / rho_safe
                vxc_contributions += direct_term
                
                # Divergence term: -(1/r²)d/dr[r²·(∂Δe/∂s)(s/|∇ρ|)sign(∇ρ)], with s/|∇ρ| = 1/(2·kf·ρ^(4/3))
                f_for_div = dE_dfeat / (2 * kf * rho_43_safe)
                div_term = cls._compute_radial_divergence(f_for_div, derivative_matrix, quadrature_nodes)
                vxc_contributions += div_term
                cls._debug_check_nan_inf(direct_term, f"direct_term for {feat_name}")
                cls._debug_check_nan_inf(div_term, f"div_term for {feat_name}")


            elif feat_name == "lap_rho_reduced":
                # Eq. 26 (appendix 3.1): q-contribution = -(5/3)(∂Δe/∂q)(q/ρ) + (1/r)d²/dr²[r·(∂Δe/∂q)(q/∇²ρ)]
                # q = ∇²ρ/(4·kf²·ρ^(5/3)) => q/∇²ρ = 1/(4·kf²·ρ^(5/3))
                rho_idx = features_list.index("rho")
                rho = features_physical[:, rho_idx]
                q = features_physical[:, feat_idx]
                kf = (3 * np.pi**2) ** (1 / 3)
                rho_safe = torch.clamp(rho, min=1e-30)
                rho_53_safe = rho_safe ** (5 / 3)
                # Direct term: -(5/3)(∂Δe/∂q)(q/ρ)
                direct_term = -(5 / 3) * dE_dfeat * q / rho_safe
                vxc_contributions += direct_term
                # Laplacian term: +(1/r)d²/dr²[r·(∂Δe/∂q)(q/∇²ρ)], with q/∇²ρ = 1/(4·kf²·ρ^(5/3))
                f_for_lap = dE_dfeat / (4 * (kf**2) * rho_53_safe)
                lap_term = cls._compute_radial_laplacian(f_for_lap, derivative_matrix, quadrature_nodes)
                vxc_contributions += lap_term
                cls._debug_check_nan_inf(direct_term, f"direct_term for {feat_name}")
                cls._debug_check_nan_inf(lap_term, f"lap_term for {feat_name}")
            
            else:
                raise ValueError(f"Feature '{feat_name}' is not currently supported for potential computation")
        
        # Debug: Check final vxc_contributions
        if cls._debug_check_nan_inf(vxc_contributions, "vxc_contributions"):
            print(f"DEBUG: NaN locations: {torch.isnan(vxc_contributions).nonzero()}")
        
        # Transform potential to transformed space if requested (for training: match y_potential_true)
        # vxc_contributions is in physical space; forward transform: symlog first, then scale
        vxc_result = vxc_contributions
        if target_component == "enhancement_factor_x":
            v_lda = cls._torch_lda_exchange_potential(rho_physical, use_safe_for_division=True)
            vxc_result = vxc_result / v_lda
        
        if transform_potential:
            # Forward transform (physical -> transformed): symlog first, then scale
            if use_symlog_for_potential and linthresh_for_potential is not None:
                # Debug: Check before symlog
                cls._debug_check_nan_inf(
                    vxc_result, "vxc_result before symlog",
                    extra=f"max(|vxc_result|) = {torch.abs(vxc_result).max()}, linthresh = {linthresh_for_potential}",
                )
                vxc_result = cls._torch_symlog(vxc_result, linthresh_for_potential)
                # Debug: Check after symlog
                cls._debug_check_nan_inf(vxc_result, "vxc_result after symlog")
            
            if scale_potential and scaler_for_potential is not None:
                # Debug: Check before scaling
                cls._debug_check_nan_inf(vxc_result, "vxc_result before scaling")
                vxc_result = cls._torch_apply_scaler(vxc_result, scaler_for_potential)
                # Debug: Check after scaling
                cls._debug_check_nan_inf(vxc_result, "vxc_result after scaling")
        
        # Reshape to match expected potential shape (n_grid,) or (n_grid, n_targets).
        # Only expand when y_energy_physical has a small, reasonable n_targets (e.g. 1 or 2 for e_x, e_c).
        # Avoid expanding when shape[1] equals n_grid (indicates model output_dim bug).
        n_grid = vxc_result.shape[0]
        if y_energy_physical.ndim > 1:
            n_targets = y_energy_physical.shape[1]
            if n_targets <= 4 and n_targets != n_grid:  # reasonable n_targets, not mistaken n_grid
                vxc_result = vxc_result.unsqueeze(-1).expand(-1, n_targets)
            else:
                vxc_result = vxc_result.unsqueeze(-1)  # (n_grid, 1) to match potential_weights
        
        return vxc_result


    def predict_potential(
        self,
        features          : np.ndarray,
        derivative_matrix : np.ndarray,
        quadrature_nodes  : np.ndarray,
        config  : ComputationConfig,
    ) -> np.ndarray:
        """
        Predict the XC potential from the features, derivative matrix, and quadrature nodes.
        """
        assert self.model_kind == "energy", "Model kind must be 'energy' for predict_potential"

        # Convert to torch tensors
        features_torch          = torch.tensor(features         , dtype=torch.float64, device=self.device)
        derivative_matrix_torch = torch.tensor(derivative_matrix, dtype=torch.float64, device=self.device)
        quadrature_nodes_torch  = torch.tensor(quadrature_nodes , dtype=torch.float64, device=self.device)

        # Compute energy density (physical) from features (physical)
        features_torch = features_torch.requires_grad_(True)
        y_energy_physical, _ = self.compute_energy_density_from_features(
            model             = self.model,
            features_physical = features_torch,
            config            = config,
            exclude_rho_from_nn = self.exclude_rho_from_nn,
        )
        
        # compute the potential
        y_potential_torch = self.compute_potential_from_energy_density(
            y_energy_physical   = y_energy_physical,
            features_physical   = features_torch,
            derivative_matrix   = derivative_matrix_torch,
            quadrature_nodes    = quadrature_nodes_torch,
            config              = config,
            requires_grad       = False,
            transform_potential = False,  # return physical potential values for inference
        )
        
        # convert the potential back to numpy; ensure 1D for single-target (n_grid,) to avoid (n_grid,1) causing downstream shape mismatch
        result = y_potential_torch.detach().cpu().numpy()
        if result.ndim == 2 and result.shape[1] == 1:
            result = result.squeeze(-1)
        return result
    

    @staticmethod
    def _compute_radial_derivative(
        f                 : torch.Tensor,
        derivative_matrix : torch.Tensor,
        quadrature_nodes  : torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute radial derivative d/dr[f] following density.py compute_density_gradient format.
        
        For spherical coordinates: d/dr[f] = [d(f·r)/dr - f] / r
        
        Parameters
        ----------
        f : torch.Tensor
            Function values at quadrature points, shape (n_grid,)
        derivative_matrix : torch.Tensor
            Derivative matrix, shape (n_elem, n_quad, n_quad) or (n_grid, n_grid)
        quadrature_nodes : torch.Tensor
            Quadrature node positions (radial coordinates), shape (n_grid,)
            
        Returns
        -------
        torch.Tensor
            Radial derivative d/dr[f], shape (n_grid,)
        """
        # Compute f·r and reshape for element-wise differentiation
        f_times_r = f * quadrature_nodes
        
        if derivative_matrix.ndim == 3:
            n_elem = derivative_matrix.shape[0]
            n_quad = derivative_matrix.shape[1]
            f_times_r_reshaped = f_times_r.reshape(n_elem, n_quad, 1)
            # Apply derivative matrix: D @ (f·r)
            d_f_r_dr = torch.bmm(derivative_matrix, f_times_r_reshaped).reshape(-1)
        else:
            # 2D derivative_matrix
            d_f_r_dr = torch.matmul(derivative_matrix, f_times_r)
        
        # Debug: Check for zero or very small quadrature_nodes
        if torch.any(torch.abs(quadrature_nodes) < 1e-10):
            zero_indices = (torch.abs(quadrature_nodes) < 1e-10).nonzero()
            print(f"DEBUG: _compute_radial_derivative: quadrature_nodes has very small values (<1e-10) at indices: {zero_indices}")
            print(f"DEBUG: min(quadrature_nodes) = {quadrature_nodes.min()}")
        
        # Compute d/dr[f] = [d(f·r)/dr - f] / r
        df_dr = (d_f_r_dr - f) / quadrature_nodes
        
        # Debug: Check for NaN/Inf in output
        TorchXCModel._debug_check_nan_inf(
            df_dr, "_compute_radial_derivative output",
            extra=f"d_f_r_dr contains NaN: {torch.isnan(d_f_r_dr).any()}, f contains NaN: {torch.isnan(f).any()}",
        )
        return df_dr


    @staticmethod
    def _compute_radial_divergence(
        f                 : torch.Tensor,
        derivative_matrix : torch.Tensor,
        quadrature_nodes  : torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute radial divergence term (1/r²) * d/dr[r² * f] for spherical coordinates.
        
        This is the divergence of a radial vector field A_r(r) r̂: ∇·(A_r r̂) = (1/r²) d/dr[r² A_r].
        Used for the GGA potential contribution from gradient terms.
        
        Parameters
        ----------
        f : torch.Tensor
            Function values at quadrature points, shape (n_grid,)
        derivative_matrix : torch.Tensor
            Derivative matrix, shape (n_elem, n_quad, n_quad) or (n_grid, n_grid)
        quadrature_nodes : torch.Tensor
            Radial coordinates, shape (n_grid,)
            
        Returns
        -------
        torch.Tensor
            (1/r²) * d/dr[r² * f], shape (n_grid,)
        """
        r_squared = quadrature_nodes ** 2
        g = r_squared * f  # r² * f
        
        if derivative_matrix.ndim == 3:
            n_elem = derivative_matrix.shape[0]
            n_quad = derivative_matrix.shape[1]
            g_reshaped = g.reshape(n_elem, n_quad, 1)
            dg_dr = torch.bmm(derivative_matrix, g_reshaped).reshape(-1)
        else:
            dg_dr = torch.matmul(derivative_matrix, g)
        
        r_squared_safe = torch.where(r_squared > 1e-20, r_squared, torch.ones_like(r_squared))
        result = dg_dr / r_squared_safe
        
        TorchXCModel._debug_check_nan_inf(result, "_compute_radial_divergence output")
        return result


    @staticmethod
    def _compute_radial_laplacian(
        f                 : torch.Tensor,
        derivative_matrix : torch.Tensor,
        quadrature_nodes  : torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute radial Laplacian nabla^2 f = (1/r)(rf)'' following density.py compute_density_laplacian.
        Here (rf)'' is the second derivative of r*f.

        Parameters
        ----------
        f : torch.Tensor
            Function values at quadrature points, shape (n_grid,)
        derivative_matrix : torch.Tensor
            Derivative matrix, shape (n_elem, n_quad, n_quad) or (n_grid, n_grid)
        quadrature_nodes : torch.Tensor
            Quadrature node positions (radial coordinates), shape (n_grid,)

        Returns
        -------
        torch.Tensor
            Radial Laplacian (1/r)(rf)'', shape (n_grid,)
        """
        if derivative_matrix.ndim == 3:
            n_elem = derivative_matrix.shape[0]
            n_quad = derivative_matrix.shape[1]

            f_reshaped = f.reshape(n_elem, n_quad, 1)
            quadrature_nodes_reshaped = quadrature_nodes.reshape(n_elem, n_quad, 1)

            # nabla^2 f = (1/r)(rf)''
            f_times_r_reshaped = quadrature_nodes_reshaped * f_reshaped
            d_rf_dr   = torch.bmm(derivative_matrix, f_times_r_reshaped)
            d2_rf_dr2 = torch.bmm(derivative_matrix, d_rf_dr)

            laplacian = (d2_rf_dr2 / quadrature_nodes_reshaped).reshape(-1)

            TorchXCModel._debug_check_nan_inf(laplacian, "_compute_radial_laplacian output")
        else:
            # 2D derivative_matrix
            f_times_r = quadrature_nodes * f
            d_rf_dr   = torch.matmul(derivative_matrix, f_times_r)
            d2_rf_dr2 = torch.matmul(derivative_matrix, d_rf_dr)
            laplacian = d2_rf_dr2 / quadrature_nodes

        return laplacian


    @staticmethod
    def _compute_weighted_inverse_laplacian(
        f                 : torch.Tensor,
        derivative_matrix : torch.Tensor,
        quadrature_nodes  : torch.Tensor,
    ) -> torch.Tensor:
        """
        Compute r * d²/dr²(f/r).

        This operator appears in the lap_rho channel of the response kernel discretization.
        Equivalent to: f'' - (2/r)f' + (2/r²)f.

        Implementation: g = f/r, dg/dr = (D@(f*r) - 2f)/r², then r*d²g/dr² = D@(dg_dr*r) - dg_dr.

        Parameters
        ----------
        f : torch.Tensor
            Function values at quadrature points, shape (n_grid,)
        derivative_matrix : torch.Tensor
            Derivative matrix, shape (n_elem, n_quad, n_quad) or (n_grid, n_grid)
        quadrature_nodes : torch.Tensor
            Quadrature node positions (radial coordinates), shape (n_grid,)

        Returns
        -------
        torch.Tensor
            r * d²/dr²(f/r), shape (n_grid,)
        """

        g = f / quadrature_nodes
        dg_dr = TorchXCModel._compute_radial_derivative(
            g, derivative_matrix, quadrature_nodes
        )
        d2g_dr2 = TorchXCModel._compute_radial_derivative(
            dg_dr, derivative_matrix, quadrature_nodes
        )
        return quadrature_nodes * d2g_dr2


    @staticmethod
    def _debug_check_nan_inf(
        x       : torch.Tensor,
        name    : str,
        context : str = "",
        extra   : str = "",
        include_stats: bool = False,
    ) -> bool:
        """
        Check tensor for NaN/Inf and print debug info if found.

        Returns True if any NaN or Inf was found, False otherwise.
        """
        has_nan = torch.isnan(x).any().item()
        has_inf = torch.isinf(x).any().item()
        if not (has_nan or has_inf):
            return False
        prefix = f"DEBUG {context}: " if context else "DEBUG: "
        if has_nan:
            print(f"{prefix}{name} contains NaN: {torch.isnan(x).sum().item()} NaNs")
        if has_inf:
            print(f"{prefix}{name} contains Inf: {torch.isinf(x).sum().item()} Infs")
        if include_stats and x.numel() > 0:
            try:
                print(f"{prefix}{name} stats: min={x.min().item()}, max={x.max().item()}, mean={x.mean().item()}")
            except Exception:
                pass
        if extra:
            print(f"{prefix}{extra}")
        return True


    @staticmethod
    def _torch_symlog(x: torch.Tensor, linthresh: float) -> torch.Tensor:
        """
        Symmetric logarithm transformation (differentiable).
        - |x| <= linthresh: y = x / linthresh
        - |x| > linthresh:  y = sign(x) * [log(|x|/linthresh) + 1]
        """
        abs_x = torch.abs(x)
        linear = abs_x <= linthresh
        out = torch.empty_like(x)
        out[linear] = x[linear] / linthresh
        out[~linear] = torch.sign(x[~linear]) * (torch.log(abs_x[~linear] / linthresh) + 1)
        return out


    @staticmethod
    def _torch_symexp(y: torch.Tensor, linthresh: float) -> torch.Tensor:
        """
        Inverse symmetric logarithm transformation (differentiable).
        - |y| <= 1: x = y * linthresh
        - |y| > 1:  x = sign(y) * linthresh * exp(|y| - 1)
        """
        abs_y = torch.abs(y)
        linear = abs_y <= 1
        out = torch.empty_like(y)
        out[linear] = y[linear] * linthresh
        exp_arg = abs_y[~linear] - 1
        exp_arg_clamped = torch.clamp(exp_arg, max=50.0)  # prevent overflow
        out[~linear] = torch.sign(y[~linear]) * linthresh * torch.exp(exp_arg_clamped)
        TorchXCModel._debug_check_nan_inf(
            out, "_torch_symexp result",
            extra=f"max(|y|)={abs_y.max()}, linthresh={linthresh}",
        )
        return out


    @staticmethod
    def _torch_symlog_prime(x: torch.Tensor, linthresh: float) -> torch.Tensor:
        """
        Derivative of symlog transformation (differentiable).
        
        Symlog'(x) = 1/linthresh if |x| <= linthresh, else sign(x)/|x|
        
        Parameters
        ----------
        x : torch.Tensor
            Input tensor (physical feature values)
        linthresh : float
            Linear threshold for symlog transformation
            
        Returns
        -------
        torch.Tensor
            Derivative of symlog at x, same shape as x
        """
        abs_x = torch.abs(x)
        symlog_prime = torch.empty_like(x)
        small_mask = abs_x <= linthresh
        symlog_prime[small_mask] = 1.0 / linthresh
        
        # Debug: Check for zero values in abs_x for large values case
        abs_x_large = abs_x[~small_mask]
        if len(abs_x_large) > 0 and torch.any(abs_x_large < 1e-12):
            zero_indices = (abs_x_large < 1e-12).nonzero()
            print(f"DEBUG: _torch_symlog_prime: abs_x has very small values (<1e-12) in large mask: {zero_indices}")
            print(f"DEBUG: min(abs_x_large) = {abs_x_large.min()}, linthresh = {linthresh}")
        
        symlog_prime[~small_mask] = torch.sign(x[~small_mask]) / abs_x[~small_mask]
        
        # Debug: Check for NaN/Inf in output
        TorchXCModel._debug_check_nan_inf(symlog_prime, "_torch_symlog_prime output")
        return symlog_prime


    @staticmethod
    def _torch_apply_scaler(x: torch.Tensor, scaler) -> torch.Tensor:
        # Get center from scaler
        if hasattr(scaler, "mean_"):
            center = scaler.mean_
        elif hasattr(scaler, "center_"):
            center = scaler.center_
        else:
            raise AttributeError(SCALER_MUST_HAVE_MEAN_OR_CENTER_ERROR.format(type(scaler).__name__))
        
        # Get scale from scaler
        if not hasattr(scaler, "scale_"):
            raise AttributeError(SCALER_MUST_HAVE_SCALE_ERROR.format(type(scaler).__name__))
        scale = scaler.scale_
        
        # Convert to numpy array first, then to torch tensor to handle numpy dtypes correctly
        center_np = np.asarray(center)
        scale_np = np.asarray(scale)


        # Ensure x.dtype is torch dtype (handle case where it might be numpy dtype)
        target_dtype = x.dtype if isinstance(x.dtype, torch.dtype) else torch.float64
        
        # Use torch.as_tensor which handles numpy arrays correctly
        # Convert to the same dtype and device as x
        center_t = torch.as_tensor(center_np, dtype=target_dtype, device=x.device)
        scale_t  = torch.as_tensor(scale_np, dtype=target_dtype, device=x.device)
        return (x - center_t) / scale_t


    @staticmethod
    def _torch_inverse_scaler(x: torch.Tensor, scaler) -> torch.Tensor:

        # Get center from scaler
        if hasattr(scaler, "mean_"):
            center = scaler.mean_
        elif hasattr(scaler, "center_"):
            center = scaler.center_
        else:
            raise AttributeError(SCALER_MUST_HAVE_MEAN_OR_CENTER_ERROR.format(type(scaler).__name__))

        # Get scale from scaler
        if not hasattr(scaler, "scale_"):
            raise AttributeError(SCALER_MUST_HAVE_SCALE_ERROR.format(type(scaler).__name__))
        scale = scaler.scale_
        
        # Convert to numpy array first, then to torch tensor to handle numpy dtypes correctly
        center_np = np.asarray(center)
        scale_np = np.asarray(scale)

        # Ensure x.dtype is torch dtype (handle case where it might be numpy dtype)
        target_dtype = x.dtype if isinstance(x.dtype, torch.dtype) else torch.float64
        
        # Use torch.as_tensor which handles numpy arrays correctly
        # Convert to the same dtype and device as x
        center_t = torch.as_tensor(center_np, dtype=target_dtype, device=x.device)
        scale_t  = torch.as_tensor(scale_np, dtype=target_dtype, device=x.device)
        return x * scale_t + center_t


    @staticmethod
    def _torch_lda_exchange_potential(
        rho: torch.Tensor,
        use_safe_for_division: bool = False,
    ) -> torch.Tensor:
        """
        LDA exchange potential v_LDA ∝ ρ^(1/3).

        Parameters
        ----------
        rho : torch.Tensor
            Electron density.
        use_safe_for_division : bool, default False
            When True, clamps the result so that |v_LDA| >= 1e-15. Use this when
            the return value will be used as a divisor (e.g., computing F_v = v_x / v_LDA)
            to avoid Inf at very small rho where v_LDA → 0.
        """
        v_lda = -0.9847450218426966 * rho ** (1 / 3)  # (3π²)^(1/3) / π
        if use_safe_for_division:
            v_lda = torch.where(
                torch.abs(v_lda) > 1e-15,
                v_lda,
                torch.full_like(v_lda, -1e-15),
            )
        return v_lda


    @staticmethod
    def _torch_lda_exchange_energy_density(
        rho: torch.Tensor,
        use_rho_43: bool = False,
    ) -> torch.Tensor:
        """
        LDA exchange energy density.

        Parameters
        ----------
        rho : torch.Tensor
            Electron density.
        use_rho_43 : bool, default False
            If True, return the full energy density directly as C * rho^(4/3).
            This avoids the gradient divergence of rho^(1/3) at small rho.
            If False, return C * rho^(1/3) (per-particle form; multiply by rho
            externally to get energy density).
        """
        if use_rho_43:
            return -0.7385587663820224 * rho ** (4 / 3)  # full epsilon_LDA_x
        return -0.7385587663820224 * rho ** (1 / 3)  # (3/4) * (3π²)^(1/3) / π


    def transform_to_torch_loader(
        self, 
        data_loader : DataLoaderType,
        batch_size  : int,
        shuffle     : bool,
        seed        : Optional[int] = None,
    ) -> torch.utils.data.DataLoader:

        # transform the data loader to a torch data loader
        from torch.utils.data import TensorDataset, DataLoader

        self.check_batch_size(batch_size)
        self.check_shuffle(shuffle)

        generator = None
        if seed is not None:
            torch.manual_seed(seed)
            generator = torch.Generator()
            generator.manual_seed(seed)

        # VxcDataLoader
        if isinstance(data_loader, VxcDataLoader):
            X = torch.DoubleTensor(data_loader.get_features_data(self.features_list))
            y = torch.DoubleTensor(data_loader.potential_transformed)
            weights = torch.DoubleTensor(data_loader.weights_for_potential)
            dataset = TensorDataset(X, y, weights)
            return DataLoader(
                dataset, 
                batch_size = batch_size, 
                shuffle    = shuffle, 
                generator  = generator
            )

        # ExcDataLoader
        elif isinstance(data_loader, ExcDataLoader):
            def collate_configs(batch):
                """
                Collate function for ExcDataLoader.
                Converts numpy arrays in ExcConfiguration instances to torch tensors,
                while preserving the list-of-dicts structure (one dict per configuration).
                Special handling for derivative_matrix: if None, use shared_derivative_matrix from data_loader.
                """
                result = []
                for cfg in batch:
                    cfg_torch = {}
                    # Convert ExcConfiguration to dict and process each field
                    cfg_dict = asdict(cfg)
                    for k, v in cfg_dict.items():
                        # Special handling for derivative_matrix
                        if k == "derivative_matrix" and v is None:
                            # Use shared_derivative_matrix from data_loader if available
                            if data_loader.shared_derivative_matrix is not None:
                                v = data_loader.shared_derivative_matrix
                            else:
                                raise ValueError("derivative_matrix is None and shared_derivative_matrix is also None.")
                        
                        if isinstance(v, np.ndarray):
                            cfg_torch[k] = torch.DoubleTensor(v)
                        else:
                            cfg_torch[k] = v  # int, None, etc. keep as-is
                    
                    # Add additional grid configuration information
                    cfg_torch["max_element_cutoff_index"] = cfg.max_element_cutoff_index
    
                    result.append(cfg_torch)
                return result  # Returns list of dict (torch tensors)
            
            # ExcDataLoader already implements __len__ and __getitem__, so can be used directly as Dataset
            return DataLoader(
                data_loader,
                batch_size = batch_size,
                shuffle    = shuffle,
                generator  = generator,
                collate_fn = collate_configs,
            )
        else:
            raise ValueError(DATA_LOADER_TYPE_NOT_VALID_ERROR.format(type(data_loader)))


    def eval_model(
        self,
        data_loader : DataLoaderType,
    ) -> Union[PotentialEvaluationMetrics, EnergyEvaluationMetrics]:
        if isinstance(data_loader, VxcDataLoader):
            return self._eval_vxc_model(data_loader)
        elif isinstance(data_loader, ExcDataLoader):
            return self._eval_exc_model(data_loader)
        else:
            raise ValueError(DATA_LOADER_TYPE_NOT_VALID_ERROR.format(type(data_loader)))



    def _eval_vxc_model(self, data_loader: VxcDataLoader) -> PotentialEvaluationMetrics:

        # check if sklearn is available
        assert SKLEARN_AVAILABLE, \
            SKLEARN_NOT_AVAILABLE_FOR_MODEL_EVALUATION_ERROR
        from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


        # check the data loader
        self._check_vxc_data_loader(data_loader)

        # get the features and target data
        self.model.eval()
        X = torch.DoubleTensor(data_loader.get_features_data(self.features_list))
        nn_indices = self._get_nn_feature_indices()
        X_for_nn = X[:, nn_indices]
        y_pred = self.model(X_for_nn).detach().numpy()
        y_true = data_loader.potential_transformed
        
        # calculate the evaluation metrics
        mae = mean_absolute_error(y_true, y_pred)
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        r2 = r2_score(y_true, y_pred)
        
        # Calculate the maximum relative error
        max_abs_true = np.max(np.abs(y_true))
        if max_abs_true > 0:
            max_relative_error = np.max(np.abs(y_true - y_pred)) / max_abs_true
        else:
            max_relative_error = np.inf if np.any(y_true != y_pred) else 0.0


        return PotentialEvaluationMetrics(
            mae  = mae,
            rmse = rmse,
            r2   = r2,
            max_relative_error = max_relative_error
        )


    def _eval_exc_model(self, data_loader: ExcDataLoader) -> EnergyEvaluationMetrics:
        """
        Evaluate energy model on ExcDataLoader.
        
        Computes metrics for energy density, potential (if include_potential=True), and total.
        """
        # check if sklearn is available
        assert SKLEARN_AVAILABLE, \
            SKLEARN_NOT_AVAILABLE_FOR_MODEL_EVALUATION_ERROR
        from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

        # check the data loader
        self._check_exc_data_loader(data_loader)
        
        include_potential        = data_loader.include_potential
        computation_config       = data_loader.get_computation_config()
        shared_derivative_matrix = data_loader.shared_derivative_matrix
        
        # Collect all predictions and true values
        self.model.eval()
        
        all_y_energy_pred = []
        all_y_energy_true = []
        all_y_potential_pred = []
        all_y_potential_true = []
        
        device = self.device if self.device else "cpu"
        for exc_configuration in data_loader.configuration_data_list:
            # Apply grid cutoff (same as _train_energy validation)
            max_element_cutoff_index = exc_configuration.max_element_cutoff_index
            quadrature_point_number  = exc_configuration.quadrature_point_number
            max_grid_cutoff_index    = max_element_cutoff_index * quadrature_point_number

            quadrature_nodes   = exc_configuration.quadrature_nodes[:max_grid_cutoff_index]
            features_physical = torch.DoubleTensor(exc_configuration.features_physical[:max_grid_cutoff_index]).to(device)
            y_energy_true     = exc_configuration.energy_density_transformed[:max_grid_cutoff_index]

            # Derivative matrix: always use shared (as in training)
            derivative_matrix = shared_derivative_matrix
            if derivative_matrix is None:
                raise ValueError("derivative_matrix is required for potential computation but shared_derivative_matrix is None")
            derivative_matrix = derivative_matrix[:max_element_cutoff_index]

            # Predict energy density via compute_energy_density_from_features
            with torch.no_grad():
                y_energy_physical, y_energy_transformed = self.compute_energy_density_from_features(
                    model             = self.model,
                    features_physical = features_physical,
                    config            = computation_config,
                    exclude_rho_from_nn = self.exclude_rho_from_nn,
                )
                y_energy_pred = y_energy_transformed.cpu().numpy()
            
            all_y_energy_pred.append(y_energy_pred)
            all_y_energy_true.append(y_energy_true)
            
            # Compute potential if include_potential is True
            if include_potential:
                quadrature_nodes_torch = torch.DoubleTensor(quadrature_nodes).to(device)
                derivative_matrix_torch = torch.DoubleTensor(derivative_matrix).to(device)
                features_with_grad = features_physical.clone().requires_grad_(True)
                
                y_energy_physical, _ = self.compute_energy_density_from_features(
                    model             = self.model,
                    features_physical = features_with_grad,
                    config            = computation_config,
                    exclude_rho_from_nn = self.exclude_rho_from_nn,
                )
                
                y_potential_pred_torch = self.compute_potential_from_energy_density(
                    y_energy_physical   = y_energy_physical,
                    features_physical   = features_with_grad,
                    derivative_matrix   = derivative_matrix_torch,
                    quadrature_nodes    = quadrature_nodes_torch,
                    config              = computation_config,
                    requires_grad       = False,
                    transform_potential = True,  # Keep in transformed space for comparison
                )
                
                y_potential_pred = y_potential_pred_torch.detach().cpu().numpy()
                y_potential_true = exc_configuration.potential_transformed[:max_grid_cutoff_index]
                
                all_y_potential_pred.append(y_potential_pred)
                all_y_potential_true.append(y_potential_true)
        
        # Concatenate all predictions and true values
        y_energy_pred_all = np.concatenate(all_y_energy_pred, axis=0)
        y_energy_true_all = np.concatenate(all_y_energy_true, axis=0)
        
        # Calculate energy metrics
        energy_mae = mean_absolute_error(y_energy_true_all.flatten(), y_energy_pred_all.flatten())
        energy_rmse = np.sqrt(mean_squared_error(y_energy_true_all.flatten(), y_energy_pred_all.flatten()))
        energy_r2 = r2_score(y_energy_true_all.flatten(), y_energy_pred_all.flatten())
        
        max_abs_energy_true = np.max(np.abs(y_energy_true_all))
        if max_abs_energy_true > 0:
            energy_max_relative_error = np.max(np.abs(y_energy_true_all - y_energy_pred_all)) / max_abs_energy_true
        else:
            energy_max_relative_error = np.inf if np.any(y_energy_true_all != y_energy_pred_all) else 0.0
        
        energy_metrics = EvaluationMetrics(
            mae                = energy_mae,
            rmse               = energy_rmse,
            r2                 = energy_r2,
            max_relative_error = energy_max_relative_error
        )
        
        # Calculate potential metrics if include_potential is True
        if include_potential:
            y_potential_pred_all = np.concatenate(all_y_potential_pred, axis=0)
            y_potential_true_all = np.concatenate(all_y_potential_true, axis=0)
            
            potential_mae = mean_absolute_error(y_potential_true_all.flatten(), y_potential_pred_all.flatten())
            potential_rmse = np.sqrt(mean_squared_error(y_potential_true_all.flatten(), y_potential_pred_all.flatten()))
            potential_r2 = r2_score(y_potential_true_all.flatten(), y_potential_pred_all.flatten())
            
            max_abs_potential_true = np.max(np.abs(y_potential_true_all))
            if max_abs_potential_true > 0:
                potential_max_relative_error = np.max(np.abs(y_potential_true_all - y_potential_pred_all)) / max_abs_potential_true
            else:
                potential_max_relative_error = np.inf if np.any(y_potential_true_all != y_potential_pred_all) else 0.0
            
            potential_metrics = EvaluationMetrics(
                mae                = potential_mae,
                rmse               = potential_rmse,
                r2                 = potential_r2,
                max_relative_error = potential_max_relative_error
            )
        else:
            potential_metrics = None
        
        # Total metrics: same as energy metrics if no potential, otherwise combine
        # For now, total = energy (as per user's requirement)
        total_metrics = energy_metrics
        
        return EnergyEvaluationMetrics(
            potential = potential_metrics,
            energy    = energy_metrics,
            total     = total_metrics
        )


    def _check_vxc_data_loader(self, data_loader: VxcDataLoader) -> None:
        # type check
        assert isinstance(data_loader, VxcDataLoader), \
            DATA_LOADER_NOT_VXCDATALOADER_ERROR.format(type(data_loader))
        for feature in self.features_list:
            assert feature in data_loader.features_list, \
                FEATURE_NOT_IN_DATA_LOADER_ERROR.format(feature, data_loader.features_list)
        if data_loader.n_samples == 0:
            raise ValueError(VXC_DATA_LOADER_CONTAINS_NO_SAMPLES_ERROR)


    def _get_nn_feature_indices(self, features_list: Optional[List[str]] = None) -> List[int]:
        """
        Return column indices of features to pass to NN.
        When exclude_rho_from_nn=True, excludes rho column; otherwise all columns.
        """
        fl = features_list if features_list is not None else self.features_list
        if self.exclude_rho_from_nn and "rho" in fl:
            return [i for i, f in enumerate(fl) if f != "rho"]
        return list(range(len(fl)))

    def _check_exc_data_loader(self, data_loader: ExcDataLoader) -> None:
        # type check
        assert isinstance(data_loader, ExcDataLoader), \
            DATA_LOADER_NOT_EXCDATALOADER_ERROR.format(type(data_loader))
        for feature in self.features_list:
            assert feature in data_loader.features_list, \
                FEATURE_NOT_IN_DATA_LOADER_ERROR.format(feature, data_loader.features_list)
        if data_loader.n_samples == 0:
            raise ValueError(EXC_DATA_LOADER_CONTAINS_NO_SAMPLES_ERROR)


    def forward(self, features: np.ndarray) -> np.ndarray:
        """
        Forward pass of the model.
        When exclude_rho_from_nn=True, only non-rho columns are passed to the NN.
        """
        # set the model to evaluation mode
        self.model.eval()

        # convert the features to a torch tensor
        features = torch.DoubleTensor(features)
        nn_indices = self._get_nn_feature_indices()
        features_for_nn = features[:, nn_indices]

        # forward pass the model, and return the predicted values
        return self.model(features_for_nn).detach().numpy()


    def _to_jsonable_model_init_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """Convert model_init_kwargs to JSON-serializable form (config scalers -> None)."""
        out = {}
        for k, v in kwargs.items():
            if isinstance(v, ComputationConfig):
                d = asdict(v)
            elif k == "config" and isinstance(v, dict):
                d = dict(v)
            else:
                out[k] = v
                continue
            for key in ("scaler_for_features", "scaler_for_potential", "scaler_for_energy"):
                if key in d:
                    d[key] = None
            out[k] = d
        return out

    def save_model(
        self,
        model_dir: Optional[Union[str, Path]] = None,
        overwrite: bool = False,
    ) -> None:  
        """
        Save the model to a directory.
        """
        # save the model
        model_dir = model_dir if model_dir is not None else self.model_dir
        model_dir = Path(model_dir)

        model_filepath = model_dir / f"{self.model_name}_model.{self.weights_ext}"
        config_filepath = model_dir / f"{self.model_name}_config.{self.config_ext}"

        if not overwrite:
            existing = [p for p in (model_filepath, config_filepath) if p.exists()]
            if existing:
                print("The following files already exist:")
                for path in existing:
                    print(f" - {path}")
                choice = input("Overwrite existing files? [y/N]: ").strip().lower()
                if choice not in ("y", "yes"):
                    print("Canceled save.")
                    return

        model_dir.mkdir(parents=True, exist_ok=True)

        torch.save(self.model.state_dict(), model_filepath)
        jsonable = self._to_jsonable_model_init_kwargs(self.model_init_kwargs)
        jsonable["exclude_rho_from_nn"] = self.exclude_rho_from_nn
        with open(config_filepath, "w") as f:
            json.dump(jsonable, f, indent=2)


    
    def compute_response_kernel_as_potential_model(
        self,
        rho                : np.ndarray,
        quadrature_nodes   : np.ndarray,
        quadrature_weights : np.ndarray,
        derivative_matrix  : np.ndarray,
        config   : ComputationConfig,
        method             : Literal["autograd", "manual"] = "autograd",
    ) -> np.ndarray:
        """
        Compute response kernel K^ρ(r_m, r_n) as a potential model via the NN Jacobian.

        Uses a discrete delta function δ_n(r) such that under the radial quadrature
        ∫₀^∞ 4πr² f(r) δ_n(r) dr ≈ f(r_n). With nodes {r_m} and weights {w_m}:
            δ_n(r_m) = δ_{mn} / (4π r_n² w_n)

        The discrete delta-perturbation ρ → ρ + h δ(·-r_n) becomes on grid samples:
            ρ_m ↦ ρ_m + h δ_n(r_m) = ρ_m + h δ_{mn} / (4π r_n² w_n)
        i.e. in vector form: ρ ↦ ρ + (h / 4π r_n² w_n) e_n

        The kernel is discretized from the NN Jacobian via the difference quotient:
            K_{mn}^ρ := K^ρ(r_m, r_n) ≈ lim_{h→0} (1/h) { v_m(ρ + (h/4πr_n²w_n)e_n) - v_m(ρ) }

        Since v_m(ρ) is smooth: v_m(ρ + α e_n) = v_m(ρ) + α ∂v_m/∂ρ_n + o(α) with α = h/(4πr_n²w_n).
        Plugging in yields the final discrete formula:
            δv(r)/δρ(r')|_{r=r_m, r'=r_n} ≈ (1 / 4π r_n² w_n) · ∂v(r_m)/∂ρ(r_n)

        Parameters
        ----------
        rho : np.ndarray
            Density on quadrature grid, shape (M,) with M = number of nodes.
        quadrature_nodes : np.ndarray
            Radial nodes {r_m}_{m=1}^M.
        quadrature_weights : np.ndarray
            Quadrature weights {w_m}_{m=1}^M.
        derivative_matrix : np.ndarray
            Unused; kept for API compatibility with energy-density-based response kernel.
        config : ComputationConfig
            Configuration for potential computation (features, scaling, symlog, etc.).
        method : {"autograd", "manual"}, default="manual"
            - "autograd": use torch.autograd.functional.jacobian on _potential_from_rho
            - "manual": use features as input and manually complete chain rule (similar to compute_potential_from_energy_density)

        Returns
        -------
        np.ndarray
            Response kernel K^ρ_{mn}, shape (M, M), where K^ρ_{mn} ≈ (1/4πr_n²w_n) ∂v_m/∂ρ_n.
        """
        if method == "autograd":
            return self.compute_response_kernel_as_potential_model_autograd(
                rho=rho,
                quadrature_nodes=quadrature_nodes,
                quadrature_weights=quadrature_weights,
                derivative_matrix=derivative_matrix,
                config=config,
            )
        else:
            return self.compute_response_kernel_as_potential_model_manual(
                rho=rho,
                quadrature_nodes=quadrature_nodes,
                quadrature_weights=quadrature_weights,
                derivative_matrix=derivative_matrix,
                config=config,
            )

    def compute_response_kernel_as_potential_model_autograd(
        self,
        rho                : np.ndarray,
        quadrature_nodes   : np.ndarray,
        quadrature_weights : np.ndarray,
        derivative_matrix  : np.ndarray,
        config   : ComputationConfig,
    ) -> np.ndarray:
        """
        Compute response kernel via autograd: jacobian(_potential_from_rho, rho).
        """
        self.model.eval()
        device = next(self.model.parameters()).device
        rho                = torch.as_tensor(rho, dtype=torch.float64, device=device)
        quadrature_nodes   = torch.as_tensor(quadrature_nodes, dtype=torch.float64, device=device)
        quadrature_weights = torch.as_tensor(quadrature_weights, dtype=torch.float64, device=device)
        derivative_matrix  = torch.as_tensor(derivative_matrix, dtype=torch.float64, device=device)

        target_component         = config.target_component_for_potential
        scale_features           = config.scale_features
        scale_potential          = config.scale_potential
        scaler_for_features      = config.scaler_for_features
        scaler_for_potential     = config.scaler_for_potential
        use_symlog_for_features  = config.use_symlog_for_features
        use_symlog_for_potential = config.use_symlog_for_potential
        linthresh_for_features   = config.linthresh_for_features
        linthresh_for_potential  = config.linthresh_for_potential

        if rho.ndim > 1:
            rho = rho.flatten()
        rho = rho.clone().detach().requires_grad_(True)

        def _potential_from_rho(rho_input: torch.Tensor) -> torch.Tensor:
            """Compute potential v(r_m) from density rho. Used for Jacobian computation."""
            features = self.compute_features_from_rho(
                rho               = rho_input,
                derivative_matrix = derivative_matrix,
                quadrature_nodes  = quadrature_nodes
            )
            features_transformed = features
            if use_symlog_for_features:
                features_transformed = self._torch_symlog(features_transformed, linthresh_for_features)
            if scale_features:
                features_transformed = self._torch_apply_scaler(features_transformed, scaler_for_features)

            nn_indices = self._get_nn_feature_indices()
            features_for_nn = features_transformed[:, nn_indices]
            potential_pred = self.model(features_for_nn)
            if potential_pred.ndim > 1:
                potential_pred = potential_pred.squeeze(-1)
            if scale_potential:
                potential_pred = self._torch_inverse_scaler(potential_pred, scaler_for_potential)
            if use_symlog_for_potential:
                potential_pred = self._torch_symexp(potential_pred, linthresh_for_potential)
            
            if target_component == "enhancement_factor_x":
                potential_pred = potential_pred * self._torch_lda_exchange_potential(rho_input)

            return potential_pred

        jacobian_v_rho = torch.autograd.functional.jacobian(
            _potential_from_rho, rho, create_graph=False
        )

        r_sq_w = 4 * np.pi * (quadrature_nodes**2) * quadrature_weights
        r_sq_w_safe = torch.clamp(r_sq_w, min=1e-20)
        scale_factor = 1.0 / r_sq_w_safe
        response_kernel = jacobian_v_rho * scale_factor.unsqueeze(1)

        return response_kernel.detach().cpu().numpy()

    def compute_response_kernel_as_potential_model_manual(
        self,
        rho                : np.ndarray,
        quadrature_nodes   : np.ndarray,
        quadrature_weights : np.ndarray,
        derivative_matrix  : np.ndarray,
        config   : ComputationConfig,
    ) -> np.ndarray:
        """
        Compute response kernel via manual chain rule: input features, manually complete
        ∂v/∂ρ = (∂v/∂features_transformed) @ (∂features_transformed/∂features) @ (∂features/∂ρ).
        Similar to compute_potential_from_energy_density.
        """
        self.model.eval()
        device = next(self.model.parameters()).device
        rho                = torch.as_tensor(rho, dtype=torch.float64, device=device)
        quadrature_nodes   = torch.as_tensor(quadrature_nodes, dtype=torch.float64, device=device)
        quadrature_weights = torch.as_tensor(quadrature_weights, dtype=torch.float64, device=device)
        derivative_matrix  = torch.as_tensor(derivative_matrix, dtype=torch.float64, device=device)

        scale_features           = config.scale_features
        scale_potential          = config.scale_potential
        scaler_for_features      = config.scaler_for_features
        scaler_for_potential     = config.scaler_for_potential
        use_symlog_for_features  = config.use_symlog_for_features
        use_symlog_for_potential = config.use_symlog_for_potential
        linthresh_for_features   = config.linthresh_for_features
        linthresh_for_potential  = config.linthresh_for_potential

        if rho.ndim > 1:
            rho = rho.flatten()
        rho = rho.clone().detach().requires_grad_(True)

        # Step 1: Compute features from rho
        features = self.compute_features_from_rho(
            rho               = rho,
            derivative_matrix = derivative_matrix,
            quadrature_nodes  = quadrature_nodes,
        )

        # Step 2: Transform features (physical → transformed)
        features_transformed = features.clone()
        if use_symlog_for_features:
            features_transformed = self._torch_symlog(features_transformed, linthresh_for_features)
        if scale_features:
            features_transformed = self._torch_apply_scaler(features_transformed, scaler_for_features)

        # Step 3: Forward pass through model
        nn_indices = self._get_nn_feature_indices()
        features_for_nn = features_transformed[:, nn_indices]
        potential_pred = self.model(features_for_nn)
        if potential_pred.ndim > 1:
            potential_pred = potential_pred.squeeze(-1)
        if scale_potential:
            potential_pred = self._torch_inverse_scaler(potential_pred, scaler_for_potential)
        if use_symlog_for_potential:
            potential_pred = self._torch_symexp(potential_pred, linthresh_for_potential)

        n_grid = rho.shape[0]
        n_features = len(self.features_list)

        # Step 4: Compute ∂v_m/∂features_transformed via autograd.grad
        dv_dfeatures_transformed = torch.autograd.grad(
            outputs      = potential_pred,
            inputs       = features_transformed,
            grad_outputs = torch.ones_like(potential_pred),
            create_graph = False,
            retain_graph = True,
            only_inputs  = True,
        )[0]

        # Step 5: Chain rule: ∂v/∂features_physical = (1/σ) * symlog' * ∂v/∂features_transformed
        features_physical = features
        if use_symlog_for_features and linthresh_for_features is not None:
            symlog_prime = self._torch_symlog_prime(
                features_physical.detach(), linthresh_for_features
            )
        else:
            symlog_prime = torch.ones_like(features_physical)
        scale_t = torch.ones(
            n_features, dtype=rho.dtype, device=rho.device
        )
        if scale_features and scaler_for_features is not None:
            scale_t = torch.as_tensor(
                scaler_for_features.scale_,
                dtype   = rho.dtype,
                device  = rho.device,
            )
        dv_dfeatures_physical = (1.0 / scale_t.unsqueeze(0)) * symlog_prime * dv_dfeatures_transformed

        # Step 6: Manually compute ∂v/∂ρ by feature type (same structure as compute_potential_from_energy_density).
        # Key: derivative operators (divergence, derivative, laplacian) act on dv_dfeat (sensitivity field),
        # NOT on ρ. v = Σ L_k(dv_dfeat_k), so ∂v/∂ρ = Σ L_k(∂dv_dfeat_k/∂ρ). Apply L to each column of ∂dv_dfeat/∂ρ.
        jacobian_v_rho = torch.zeros(
            n_grid, n_grid, dtype=rho.dtype, device=rho.device
        )

        for feat_idx, feat_name in enumerate(self.features_list):
            dv_dfeat = dv_dfeatures_physical[:, feat_idx]  # (n_grid,) = ∂v/∂feat

            if feat_name == "rho":
                # Direct: v += dv_dfeat, so ∂v/∂ρ += diag(dv_dfeat)
                jacobian_v_rho += torch.diag(dv_dfeat)

            elif feat_name in ["grad_rho_norm", "grad_rho_abs", "grad_rho_mag"]:
                # v += -sign(∇ρ) * d/dr[dv_dfeat] = + _compute_radial_derivative(dv_dfeat)
                for n in range(n_grid):
                    jacobian_v_rho[:, n] += self._compute_radial_derivative(
                        dv_dfeat, derivative_matrix, quadrature_nodes
                    )

            elif feat_name == "lap_rho":
                # v += laplacian(dv_dfeat), so ∂v/∂ρ = L_lap(∂dv_dfeat/∂ρ)
                for n in range(n_grid):
                    jacobian_v_rho[:, n] += self._compute_weighted_inverse_laplacian(
                        dv_dfeat, derivative_matrix, quadrature_nodes
                    )

            elif feat_name in ["grad_rho_reduced", "lap_rho_reduced"]:
                raise NotImplementedError(
                    f"Feature '{feat_name}' is not yet supported in "
                    "compute_response_kernel_as_potential_model_manual. "
                    "Supported: rho, grad_rho, grad_rho_norm, grad_rho_abs, grad_rho_mag, lap_rho."
                )
            else:
                raise ValueError(INVALID_FEATURE_NAME_FOR_RESPONSE_KERNEL_COMPUTATION_ERROR.format(feat_name))

        r_sq_w = 4 * np.pi * (quadrature_nodes**2) * quadrature_weights
        r_sq_w_safe = torch.clamp(r_sq_w, min=1e-20)
        scale_factor = 1.0 / r_sq_w_safe
        response_kernel = jacobian_v_rho * scale_factor.unsqueeze(1)

        return response_kernel.detach().cpu().numpy()





    def compute_features_from_rho(
        self, 
        rho               : torch.Tensor,
        derivative_matrix : torch.Tensor,
        quadrature_nodes  : torch.Tensor,
    ) -> torch.Tensor:

        features = []
        for feature_name in self.features_list:
            if feature_name == "rho":
                features.append(rho)
            elif feature_name == "grad_rho":
                grad_rho = self._compute_radial_derivative(rho, derivative_matrix, quadrature_nodes)
                features.append(grad_rho)
            elif feature_name in ["grad_rho_norm", "grad_rho_abs", "grad_rho_mag"]:
                grad_rho = self._compute_radial_derivative(rho, derivative_matrix, quadrature_nodes)
                features.append(torch.abs(grad_rho))
            elif feature_name == "lap_rho":
                lap_rho = self._compute_radial_laplacian(rho, derivative_matrix, quadrature_nodes)
                features.append(lap_rho)
            elif feature_name == "grad_rho_reduced":
                grad_rho = self._compute_radial_derivative(rho, derivative_matrix, quadrature_nodes)
                kf = (3 * np.pi**2) ** (1 / 3)
                grad_rho_reduced = torch.abs(grad_rho) / (2 * kf * rho ** (4 / 3))
                features.append(grad_rho_reduced)
            elif feature_name == "lap_rho_reduced":
                lap_rho = self._compute_radial_laplacian(rho, derivative_matrix, quadrature_nodes)
                kf = (3 * np.pi**2) ** (1 / 3)
                lap_rho_reduced = lap_rho / (4 * (kf**2) * rho ** (5 / 3))
                features.append(lap_rho_reduced)
            else:
                raise ValueError(INVALID_FEATURE_NAME_FOR_RESPONSE_KERNEL_COMPUTATION_ERROR.format(feature_name))
        
        features = torch.stack(features, dim=1)
        return features



