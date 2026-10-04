"""
Calculator class to hold a model plus its pipeline metadata.
"""

from __future__ import annotations


import json
import pickle
import numpy as np
import importlib

from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union

from .base import XCBaseModel
from ...data.data_processing import DataProcessor
from ...data.data_manager import (
    VxcDataLoader,
    ExcDataLoader,
    ComputationConfig,
    LossNorm,
)

DataLoaderType = Union[VxcDataLoader, ExcDataLoader]

try:
    from sklearn.preprocessing import StandardScaler, RobustScaler
    SKLEARN_AVAILABLE = True
    ScalerType = Union[StandardScaler, RobustScaler]
except ImportError:
    SKLEARN_AVAILABLE = False
    ScalerType = Any


# Error messages
LINTHRESH_TARGETS_NOT_FOUND_ERROR = \
    "linthresh_for_potential is required when use_symlog_for_potential is True."
LINTHRESH_ENERGY_NOT_FOUND_ERROR = \
    "linthresh_for_energy is required when use_symlog_for_energy is True."
LINTHRESH_FEATURES_NOT_FOUND_ERROR = \
    "linthresh_for_features is required when use_symlog_for_features is True."

CALLING_PREDICT_EXC_FOR_MODEL_KIND_POTENTIAL_ERROR = \
    "Model kind is potential, but predict_exc is called."
CALLING_COMPUTE_ENERGY_FOR_MODEL_KIND_POTENTIAL_ERROR = \
    "Model kind is potential, but compute_energy is called."

DERIVATIVE_MATRIX_AND_QUADRATURE_NODES_REQUIRED_ERROR = \
    "derivative_matrix and quadrature_nodes are required when model_kind is 'energy'."


@dataclass
class MLXCCalculator:
    """
    Backend-agnostic wrapper for model + preprocessing metadata.
    """

    model                       : XCBaseModel
    features_list               : List[str]
    target_functional           : str
    target_component            : str
    target_mode                 : str
    reference_functional        : Optional[str]
    
    # required parameters for scaling (no defaults)
    scale_features              : bool
    scale_potential             : bool
    scale_energy                : Optional[bool]
    scaler_type_for_features    : str
    scaler_type_for_potential   : str
    scaler_type_for_energy      : Optional[str]

    # optional parameters (with defaults)
    target_component_for_potential : Optional[str] = None  # When set (e.g. "v_xc"), used for potential config instead of target_component
    scaler_kwargs_for_features  : Dict[str, Any] = field(default_factory=dict)
    scaler_kwargs_for_potential : Dict[str, Any] = field(default_factory=dict)
    scaler_kwargs_for_energy    : Dict[str, Any] = field(default_factory=dict)
    scaler_for_features         : Optional[ScalerType] = None
    scaler_for_potential        : Optional[ScalerType] = None
    scaler_for_energy           : Optional[ScalerType] = None
    
    # optional parameters for symlog transformation
    use_symlog_for_features     : bool = True
    use_symlog_for_potential    : bool = True
    use_symlog_for_energy       : Optional[bool]  = None
    linthresh_for_features      : Optional[float] = None
    linthresh_for_potential     : Optional[float] = None
    linthresh_for_energy        : Optional[float] = None

    # exclude rho from NN input (for enhancement factor with reduced features only)
    exclude_rho_from_nn         : bool = False

    # sample-weight convention from dataloaders (L1norm=MAE weights, L2norm=MSE weights)
    loss_norm_for_potential     : LossNorm = "L1norm"
    loss_norm_for_energy        : Optional[LossNorm] = None



    def print_info(self, print_footer_separator: bool = True) -> None:
        """
        Print information about the calculator.
        """
        print("=" * 75)
        print("Machine Learning XC Calculator".center(75))
        print("=" * 75)
        print(f"\t model                     : {self.model.model_name}")
        print(f"\t features_list             : {self.features_list}")
        print(f"\t target_functional         : {self.target_functional}")
        print(f"\t target_component          : {self.target_component}")
        print(f"\t target_mode               : {self.target_mode}")
        print(f"\t reference_functional      : {self.reference_functional}")
        # parameters for scaling
        print(f"\t scale_features            : {self.scale_features}")
        print(f"\t scale_potential           : {self.scale_potential}")
        print(f"\t scale_energy              : {self.scale_energy}")
        print(f"\t scaler_type_for_features  : {self.scaler_type_for_features}")
        print(f"\t scaler_type_for_potential : {self.scaler_type_for_potential}")
        print(f"\t scaler_type_for_energy    : {self.scaler_type_for_energy}")
        print(f"\t scaler_for_features       : {self.scaler_for_features}")
        print(f"\t scaler_for_potential      : {self.scaler_for_potential}")
        print(f"\t scaler_for_energy         : {self.scaler_for_energy}")
        # parameters for symlog transformation
        print(f"\t use_symlog_for_features   : {self.use_symlog_for_features}")
        print(f"\t use_symlog_for_potential  : {self.use_symlog_for_potential}")
        print(f"\t use_symlog_for_energy     : {self.use_symlog_for_energy}")
        print(f"\t linthresh_for_features    : {self.linthresh_for_features}")
        print(f"\t linthresh_for_potential   : {self.linthresh_for_potential}")
        print(f"\t linthresh_for_energy      : {self.linthresh_for_energy}")
        print(f"\t exclude_rho_from_nn       : {self.exclude_rho_from_nn}")
        print(f"\t loss_norm_for_potential    : {self.loss_norm_for_potential}")
        if self.loss_norm_for_energy is not None:
            print(f"\t loss_norm_for_energy      : {self.loss_norm_for_energy}")
        print()
        
        if print_footer_separator:
            print("=" * 75)



    @classmethod
    def from_dataloader(
        cls,
        model       : XCBaseModel,
        data_loader : DataLoaderType,
    ) -> "MLXCCalculator":
        # Check type using isinstance first, then fallback to type name check
        # This handles cases where imports might come from different paths
        if isinstance(data_loader, VxcDataLoader):
            return cls.from_vxc_dataloader(model, data_loader)
        elif isinstance(data_loader, ExcDataLoader):
            return cls.from_exc_dataloader(model, data_loader)
        else:
            # Fallback: check by type name to handle import path differences
            data_loader_type_name = type(data_loader).__name__
            if data_loader_type_name == "VxcDataLoader":
                return cls.from_vxc_dataloader(model, data_loader)
            elif data_loader_type_name == "ExcDataLoader":
                return cls.from_exc_dataloader(model, data_loader)
            else:
                raise ValueError(f"Unsupported data loader type: {type(data_loader)}")


    @classmethod
    def from_vxc_dataloader(
        cls,
        model       : XCBaseModel,
        data_loader : DataLoaderType,
    ) -> "MLXCCalculator":

        return cls(
            model                       = model,
            features_list               = data_loader.features_list,
            target_functional           = data_loader.target_functional,
            target_component            = data_loader.target_component,
            target_mode                 = data_loader.target_mode,
            reference_functional        = data_loader.reference_functional,
            # parameters for scaling
            scale_features              = data_loader.scale_features,
            scale_potential             = data_loader.scale_potential,
            scale_energy                = None,
            scaler_type_for_features    = data_loader.scaler_type_for_features,
            scaler_type_for_potential   = data_loader.scaler_type_for_potential,
            scaler_type_for_energy      = None,
            scaler_kwargs_for_features  = data_loader.scaler_kwargs_for_features,
            scaler_kwargs_for_potential = data_loader.scaler_kwargs_for_potential,
            scaler_kwargs_for_energy    = {},
            scaler_for_features         = data_loader.scaler_for_features,
            scaler_for_potential        = data_loader.scaler_for_potential,
            scaler_for_energy           = None,
            # parameters for symlog transformation
            use_symlog_for_features     = data_loader.use_symlog_for_features,
            use_symlog_for_potential    = data_loader.use_symlog_for_potential,
            use_symlog_for_energy       = None,
            linthresh_for_features      = data_loader.linthresh_for_features,
            linthresh_for_potential     = data_loader.linthresh_for_potential,
            linthresh_for_energy        = None,
            exclude_rho_from_nn         = getattr(data_loader, "exclude_rho_from_nn", getattr(model, "exclude_rho_from_nn", False)),
            loss_norm_for_potential     = data_loader.loss_norm_for_potential,
            loss_norm_for_energy        = None,
        )


    @classmethod
    def from_exc_dataloader(
        cls,
        model       : XCBaseModel,
        data_loader : DataLoaderType,
    ) -> "MLXCCalculator":

        return cls(
            model                       = model,    
            features_list               = data_loader.features_list,
            target_functional           = data_loader.target_functional,
            target_component            = data_loader.target_component,
            target_component_for_potential = getattr(data_loader, "target_component_for_potential", None),
            target_mode                 = data_loader.target_mode,
            reference_functional        = data_loader.reference_functional,
            # parameters for scaling
            scale_features              = data_loader.scale_features,
            scale_potential             = data_loader.scale_potential,
            scale_energy                = data_loader.scale_energy,
            scaler_type_for_features    = data_loader.scaler_type_for_features,
            scaler_type_for_potential   = data_loader.scaler_type_for_potential,
            scaler_type_for_energy      = data_loader.scaler_type_for_energy,
            scaler_kwargs_for_features  = data_loader.scaler_kwargs_for_features,
            scaler_kwargs_for_potential = data_loader.scaler_kwargs_for_potential,
            scaler_kwargs_for_energy    = data_loader.scaler_kwargs_for_energy,
            scaler_for_features         = data_loader.scaler_for_features,
            scaler_for_potential        = data_loader.scaler_for_potential,
            scaler_for_energy           = data_loader.scaler_for_energy,
            # parameters for symlog transformation
            use_symlog_for_features     = data_loader.use_symlog_for_features,
            use_symlog_for_potential    = data_loader.use_symlog_for_potential,
            use_symlog_for_energy       = data_loader.use_symlog_for_energy,
            linthresh_for_features      = data_loader.linthresh_for_features,
            linthresh_for_potential     = data_loader.linthresh_for_potential,
            linthresh_for_energy        = data_loader.linthresh_for_energy,
            exclude_rho_from_nn         = getattr(data_loader, "exclude_rho_from_nn", getattr(model, "exclude_rho_from_nn", False)),
            loss_norm_for_potential     = data_loader.loss_norm_for_potential,
            loss_norm_for_energy        = data_loader.loss_norm_for_energy,
        )



    @classmethod
    def load(
        cls,
        model_dir    : Union[str, Path],
        model_name   : str,
        model_cls    : Optional[type] = None,
        model_kind   : Optional[str] = None,
        features_list: Optional[List[str]] = None,
        device       : Optional[str] = "cpu",
    ) -> "MLXCCalculator":
        model_dir = Path(model_dir)
        config_path = model_dir / f"{model_name}_config.json"
        scalers_path = model_dir / f"{model_name}_scalers.pkl"
        metadata_path = model_dir / f"{model_name}_metadata.json"

        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        if not scalers_path.exists():
            raise FileNotFoundError(f"Scalers file not found: {scalers_path}")
        if not metadata_path.exists():
            raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

        with config_path.open("r", encoding="utf-8") as f:
            model_config = json.load(f)
        if isinstance(model_config, dict) and "model_config" in model_config: # This is a fix for the old model config format
            model_config = model_config["model_config"]
        with scalers_path.open("rb") as f:
            scalers = pickle.load(f)
        with metadata_path.open("r", encoding="utf-8") as f:
            metadata = json.load(f)

        if model_cls is None:
            model_class_path = metadata.get("model_class")
            if model_class_path is None:
                raise ValueError("'model_cls' is required when metadata does not include model_class.")
            module_path, class_name = model_class_path.rsplit(".", 1)
            module = importlib.import_module(module_path)
            model_cls = getattr(module, class_name)

        if model_kind is None:
            model_kind = metadata.get("model_kind")
        if features_list is None:
            features_list = metadata.get("features_list", [])

        # Extract exclude_rho_from_nn from config or metadata (default False for backward compatibility)
        exclude_rho_from_nn = model_config.pop("exclude_rho_from_nn", metadata.get("exclude_rho_from_nn", False))

        from .torch_backend import TorchXCModel
        model = TorchXCModel(
            model_kind        = model_kind,
            features_list     = features_list,
            model_cls         = model_cls,
            model_init_kwargs = model_config,
            model_name        = model_name,
            model_dir         = str(model_dir),
            device            = device,
            exclude_rho_from_nn = exclude_rho_from_nn,
        )

        model_file = model_dir / f"{model_name}_model.{model.weights_ext}"
        if not model_file.exists():
            raise FileNotFoundError(f"Model file not found: {model_file}")
        import torch
        state_dict = torch.load(model_file, map_location=device)
        model.model.load_state_dict(state_dict)

        # Inject scalers into model's config dict (for energy models using config.to_dict())
        if hasattr(model.model, "config") and isinstance(model.model.config, dict):
            model.model.config["scaler_for_features"] = scalers.get("scaler_for_features") or scalers.get("scaler_X")
            model.model.config["scaler_for_potential"] = scalers.get("scaler_for_potential") or scalers.get("scaler_y")
            model.model.config["scaler_for_energy"] = scalers.get("scaler_for_energy")

        return cls(
            model                       = model,
            features_list               = metadata.get("features_list", []),
            target_functional           = metadata.get("target_functional", ""),
            target_component            = metadata.get("target_component", ""),
            target_component_for_potential = metadata.get("target_component_for_potential", None),
            target_mode                 = metadata.get("target_mode", ""),
            reference_functional        = metadata.get("reference_functional", None),
            scale_features              = metadata.get("scale_features", True),
            scale_potential             = metadata.get("scale_potential", metadata.get("scale_targets", True)),
            scale_energy                = metadata.get("scale_energy", None),
            scaler_type_for_features    = metadata.get("scaler_type_for_features",  metadata.get("scaler_type_features", "robust")),
            scaler_type_for_potential   = metadata.get("scaler_type_for_potential", metadata.get("scaler_type_targets", "robust")),
            scaler_type_for_energy      = metadata.get("scaler_type_for_energy", None),
            scaler_kwargs_for_features  = metadata.get("scaler_kwargs_for_features",  metadata.get("scaler_kwargs_features", {})),
            scaler_kwargs_for_potential = metadata.get("scaler_kwargs_for_potential", metadata.get("scaler_kwargs_targets", {})),
            scaler_kwargs_for_energy    = metadata.get("scaler_kwargs_for_energy", {}),
            scaler_for_features         = scalers.get("scaler_for_features")  if "scaler_for_features"  in scalers else scalers.get("scaler_X"),
            scaler_for_potential        = scalers.get("scaler_for_potential") if "scaler_for_potential" in scalers else scalers.get("scaler_y"),
            scaler_for_energy           = scalers.get("scaler_for_energy", None),
            use_symlog_for_features     = metadata.get("use_symlog_for_features",  metadata.get("use_symlog_features", True)),
            use_symlog_for_potential    = metadata.get("use_symlog_for_potential", metadata.get("use_symlog_targets",  True)),
            use_symlog_for_energy       = metadata.get("use_symlog_for_energy", None),
            linthresh_for_features      = metadata.get("linthresh_for_features",  metadata.get("linthresh_features", None)),
            linthresh_for_potential     = metadata.get("linthresh_for_potential", metadata.get("linthresh_targets",  None)),
            linthresh_for_energy        = metadata.get("linthresh_for_energy", None),
            exclude_rho_from_nn         = metadata.get("exclude_rho_from_nn", False),
            loss_norm_for_potential     = metadata.get("loss_norm_for_potential", "L1norm"),
            loss_norm_for_energy        = metadata.get("loss_norm_for_energy", None),
        )


    def save(
        self, 
        model_dir: Optional[Union[str, Path]] = None,
        overwrite: bool = False
    ) -> None:

        """
        Save the calculator to a directory.
        
        Parameters
        ----------
        model_dir : Optional[Union[str, Path]]
            Directory to save the calculator.
        """

        model_dir = model_dir if model_dir is not None else self.model.model_dir
        model_dir = Path(model_dir)

        if not overwrite:
            existing = [p for p in self._resolve_model_paths(model_dir) if p.exists()]
            if existing:
                print("Saving MLXC calculator detected existing files.")
                print("The following files already exist:")
                for path in existing:
                    print(f" - {path}")
                choice = input("Overwrite existing files? [y/N]: ").strip().lower()
                if choice not in ("y", "yes"):
                    print("Canceled save.")
                    return

        model_dir.mkdir(parents=True, exist_ok=True)

        # save the model
        self.model.save_model(model_dir, overwrite=True)


        scalers_path = self._resolve_scalers_path(model_dir)
        scalers_data = {
            "scaler_X": self.scaler_for_features,
            "scaler_y": self.scaler_for_potential,
            "scaler_for_energy": self.scaler_for_energy,
        }
        with scalers_path.open("wb") as f:
            pickle.dump(scalers_data, f)

        metadata_path = self._resolve_metadata_path(model_dir)
        metadata_json = self._jsonable(self._to_metadata())
        with metadata_path.open("w", encoding="utf-8") as f:
            json.dump(metadata_json, f, indent=4)



    def inverse_transform_potential_predictions(self, y_pred: np.ndarray) -> np.ndarray:
        """
        Inverse transform potential predictions back to physical space.
        """
        y_pred = np.asarray(y_pred).copy()

        # Step 1: Inverse scaling (if scaler was used)
        if self.scaler_for_potential is not None:
            y_pred = self.scaler_for_potential.inverse_transform(y_pred.reshape(-1, 1)).flatten()
        
        # Step 2: Inverse symlog (if symlog was used)
        if self.use_symlog_for_potential:
            assert self.linthresh_for_potential is not None, \
                LINTHRESH_TARGETS_NOT_FOUND_ERROR
            y_pred = DataProcessor.symexp(y_pred, linthresh=self.linthresh_for_potential)

        return y_pred


    def inverse_transform_energy_predictions(self, y_pred: np.ndarray) -> np.ndarray:
        """
        Inverse transform energy predictions back to physical space.
        Uses scaler_for_energy and use_symlog_for_energy (not potential).
        """
        y_pred = np.asarray(y_pred).copy()

        # Step 1: Inverse scaling (if scaler was used)
        if self.scaler_for_energy is not None:
            y_pred = self.scaler_for_energy.inverse_transform(y_pred.reshape(-1, 1)).flatten()

        # Step 2: Inverse symlog (if symlog was used)
        if self.use_symlog_for_energy:
            assert self.linthresh_for_energy is not None, LINTHRESH_ENERGY_NOT_FOUND_ERROR
            y_pred = DataProcessor.symexp(y_pred, linthresh=self.linthresh_for_energy)

        return y_pred


    def transform_features(self, X: np.ndarray) -> np.ndarray:
        """
        Transform features to the transformed space.
        """
        X = np.asarray(X).copy()
        

        # Step 1: Transform features (if symlog was used)
        if self.use_symlog_for_features:
            assert self.linthresh_for_features is not None, \
                LINTHRESH_FEATURES_NOT_FOUND_ERROR
            X = DataProcessor.symlog(X, linthresh=self.linthresh_for_features)

        # Step 2: Transform features (if scaler was used)
        if self.scaler_for_features is not None:
            X = self.scaler_for_features.transform(X)
        

        return X


    def inverse_transform_features(self, X_transformed: np.ndarray) -> np.ndarray:
        """
        Inverse transform features back to physical space.
        """
        X_transformed = np.asarray(X_transformed).copy()
        
        # Step 1: Inverse scaling (if scaler was used)
        if self.scaler_for_features is not None:
            X_transformed = self.scaler_for_features.inverse_transform(X_transformed)
        
        # Step 2: Inverse symlog (if symlog was used)
        if self.use_symlog_for_features:
            assert self.linthresh_for_features is not None, \
                LINTHRESH_FEATURES_NOT_FOUND_ERROR
            X_transformed = DataProcessor.symexp(X_transformed, linthresh=self.linthresh_for_features)
        
        return X_transformed

    
    def predict_vxc(
        self, 
        features          : np.ndarray,
        derivative_matrix : Optional[np.ndarray] = None,
        quadrature_nodes  : Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Forward pass of the calculator.
        
        Parameters
        ----------
        features : np.ndarray
            Input features, shape (n_grid, n_features).
        derivative_matrix : Optional[np.ndarray], default=None
            Derivative matrix for computing gradients in spherical coordinates.
            Required when model_kind == "energy".
            Shape: (n_elements, n_quad_per_element, n_quad_per_element) or (n_grid, n_grid).
        quadrature_nodes : Optional[np.ndarray], default=None
            Quadrature node positions (radial coordinates), shape (n_grid,).
            Required when model_kind == "energy".
        
        Returns
        -------
        np.ndarray
            Predicted XC potential values, shape (n_grid,) or (n_grid, n_targets).
        """
        if self.model_kind == "potential":
            features = np.asarray(features).copy()
            features_transformed = self.transform_features(features)
            predictions = self.model.forward(features_transformed)
            predictions = self.inverse_transform_potential_predictions(predictions)
            
            # If target component is enhancement factor x, convert F_v to v_x = F_v * v_LDA
            if self.target_component == "enhancement_factor_x":
                from ..lda import lda_exchange_potential_generic
                rho = features[:, self.features_list.index("rho")]
                lda_exchange_potential = lda_exchange_potential_generic(rho)
                # Match lda shape to predictions for element-wise multiply (avoid (n,1)*(n,) -> (n,n) broadcast)
                if np.asarray(predictions).ndim > 1:
                    lda_exchange_potential = lda_exchange_potential.reshape(-1, 1)
                predictions = predictions * lda_exchange_potential
                
            return predictions
        else:
            # Model kind is "energy": predict energy density first, then compute potential
            if derivative_matrix is None or quadrature_nodes is None:
                raise ValueError(DERIVATIVE_MATRIX_AND_QUADRATURE_NODES_REQUIRED_ERROR)
            
            # predict the potential (returns raw v_xc when potential target is v_xc, or F_v when enhancement_factor_x)
            computation_config = self.get_computation_config()
            predictions = self.model.predict_potential(
                features          = features,
                derivative_matrix = derivative_matrix,
                quadrature_nodes  = quadrature_nodes,
                config            = computation_config,
            )
            
            return predictions

    def predict_vxc_with_spin(self, features: np.ndarray) -> np.ndarray:
        raise NotImplementedError("predict_vxc_with_spin is not implemented for ML XC models")


    def predict_exc(self, features: np.ndarray) -> np.ndarray:
        """
        Forward pass of the calculator.
        Uses energy scaler and symlog (inverse_transform_energy_predictions),
        not the potential scaler.

        For target_component == "enhancement_factor_x", the model predicts F (enhancement factor).
        Physical energy density is e = F * (ε_LDA * ρ). This method returns e in physical space.
        """
        if self.model_kind == "potential":
            raise ValueError(CALLING_PREDICT_EXC_FOR_MODEL_KIND_POTENTIAL_ERROR)

        # Extract rho in physical space before transforming (needed for enhancement_factor_x)
        rho_physical = features[:, self.features_list.index("rho")] if self.target_component == "enhancement_factor_x" else None

        features    = self.transform_features(features)
        predictions = self.model.forward(features)
        predictions = self.inverse_transform_energy_predictions(predictions)

        # If target is enhancement factor, convert F -> physical energy density: e = F * (ε_LDA * ρ)
        if self.target_component == "enhancement_factor_x":
            from ..lda import lda_exchange_energy_density_generic
            lda_energy_density = lda_exchange_energy_density_generic(rho_physical) * rho_physical
            if np.asarray(predictions).ndim > 1:
                lda_energy_density = lda_energy_density.reshape(-1, 1)
            predictions = predictions * lda_energy_density

        return predictions



    def compute_response_kernel(
        self, 
        rho                : np.ndarray,
        quadrature_nodes   : np.ndarray,
        quadrature_weights : np.ndarray,
        derivative_matrix  : np.ndarray,
    ) -> np.ndarray:
        """
        Response kernel.
        """
        if self.model_kind == "potential":
            return self.model.compute_response_kernel_as_potential_model(
                rho                = rho,
                quadrature_nodes   = quadrature_nodes,
                quadrature_weights = quadrature_weights,
                derivative_matrix  = derivative_matrix,
                config             = self.get_computation_config(),
            )
        else:
            raise NotImplementedError("Response kernel is not implemented for ML XC models")



    def compute_energy(
        self,
        density_data       : Any,
        quadrature_nodes   : np.ndarray,
        quadrature_weights : np.ndarray,
        density_calculator : Optional[Any] = None,
        poisson_solver     : Optional[Any] = None,
    ) -> float:
        """Compute ML XC energy: construct features (driver-style) -> predict_exc -> integrate."""
        if self.model.model_kind == "potential":
            raise ValueError(CALLING_COMPUTE_ENERGY_FOR_MODEL_KIND_POTENTIAL_ERROR)

        features = self._construct_features(density_data, density_calculator, poisson_solver)
        e_density = self.predict_exc(features).flatten()
        return float(np.sum(4 * np.pi * quadrature_nodes**2 * quadrature_weights * e_density))


    def _construct_features(
        self,
        density_data       : Any,
        density_calculator : Optional[Any] = None,
        poisson_solver     : Optional[Any] = None,
    ) -> np.ndarray:
        """Construct features for ML model (same logic as SCFDriver._construct_input_features_for_ml_model)."""
        rho = np.asarray(density_data.rho).flatten()
        features = []
        for feat in self.features_list:
            if feat == "rho":
                features.append(rho)
            elif feat == "grad_rho":
                grad_rho = density_data.grad_rho if density_data.grad_rho is not None else (density_calculator.compute_density_gradient(rho) if density_calculator else np.zeros_like(rho))
                features.append(np.asarray(grad_rho).flatten())
            elif feat == "grad_rho_norm":
                grad_rho = density_data.grad_rho if density_data.grad_rho is not None else (density_calculator.compute_density_gradient(rho) if density_calculator else np.zeros_like(rho))
                features.append(np.abs(np.asarray(grad_rho).flatten()))
            elif feat == "lap_rho":
                lap_rho = density_calculator.compute_density_laplacian(rho) if density_calculator else np.zeros_like(rho)
                features.append(lap_rho)
            elif feat == "grad_rho_reduced":
                grad_rho = density_data.grad_rho if density_data.grad_rho is not None else (density_calculator.compute_density_gradient(rho) if density_calculator else np.zeros_like(rho))
                kf = (3 * np.pi**2) ** (1.0 / 3.0)
                rho_safe = np.maximum(rho, 1e-20)
                features.append(np.abs(np.asarray(grad_rho).flatten()) / (2.0 * kf * (rho_safe ** (4.0 / 3.0))))
            elif feat == "lap_rho_reduced":
                lap_rho = density_calculator.compute_density_laplacian(rho) if density_calculator else np.zeros_like(rho)
                kf = (3 * np.pi**2) ** (1.0 / 3.0)
                rho_safe = np.maximum(rho, 1e-20)
                features.append(lap_rho / (4.0 * (kf ** 2) * (rho_safe ** (5.0 / 3.0))))
            elif feat == "hartree":
                hartree = poisson_solver.solve_hartree(rho) if poisson_solver else np.zeros_like(rho)
                features.append(hartree)
            elif feat == "lda_xc_potential":
                from ..lda import lda_exchange_potential_generic, lda_correlation_potential_generic
                lda_xc = lda_exchange_potential_generic(rho) + lda_correlation_potential_generic(rho)
                features.append(lda_xc)
            else:
                raise ValueError(f"compute_energy: unsupported feature '{feat}'")
        return np.column_stack(features)


    def get_computation_config(self) -> ComputationConfig:
        """
        Get configuration for energy and potential computation.
        
        Returns a ComputationConfig dataclass containing all necessary parameters
        for computing energy density and XC potential from model predictions.
        
        Returns
        -------
        ComputationConfig
            Configuration object with all parameters needed for computation.
        """
        target_comp_pot = self.target_component_for_potential if self.target_component_for_potential is not None else self.target_component

        return ComputationConfig(
            features_list                  = self.features_list,
            target_component_for_potential = target_comp_pot,
            target_component_for_energy    = self.target_component,
            # optional parameters for scaling
            scale_features           = self.scale_features,
            scale_potential          = self.scale_potential,
            scale_energy             = self.scale_energy,
            scaler_for_features      = self.scaler_for_features,
            scaler_for_potential     = self.scaler_for_potential,
            scaler_for_energy        = self.scaler_for_energy,
            # optional parameters for symlog transformation
            use_symlog_for_features  = self.use_symlog_for_features,
            use_symlog_for_potential = self.use_symlog_for_potential,
            use_symlog_for_energy    = self.use_symlog_for_energy,
            linthresh_for_features   = self.linthresh_for_features,
            linthresh_for_potential  = self.linthresh_for_potential,
            linthresh_for_energy     = self.linthresh_for_energy,
        )

    # Backward compatibility alias
    get_potential_computation_config = get_computation_config

    def _to_metadata(self) -> Dict[str, Any]:
        model_class = None
        if getattr(self.model, "model_cls", None) is not None:
            model_class = f"{self.model.model_cls.__module__}.{self.model.model_cls.__name__}"
        return {
            # parameters for documentation
            "features_list"                : self.features_list,
            "target_functional"            : self.target_functional,
            "target_component"             : self.target_component,
            "target_component_for_potential": self.target_component_for_potential,
            "target_mode"                  : self.target_mode,
            "reference_functional"         : self.reference_functional,
            "model_kind"                   : getattr(self.model, "model_kind", None),
            "model_class"                  : model_class,
            
            # parameters for scaling
            "scale_features"               : self.scale_features,
            "scale_potential"              : self.scale_potential,
            "scale_energy"                 : self.scale_energy,
            "scaler_type_for_features"     : self.scaler_type_for_features,
            "scaler_type_for_potential"    : self.scaler_type_for_potential,
            "scaler_type_for_energy"       : self.scaler_type_for_energy,
            "scaler_kwargs_for_features"   : self.scaler_kwargs_for_features,
            "scaler_kwargs_for_potential"  : self.scaler_kwargs_for_potential,
            "scaler_kwargs_for_energy"     : self.scaler_kwargs_for_energy,
            
            # parameters for symlog transformation
            "use_symlog_for_features"      : self.use_symlog_for_features,
            "use_symlog_for_potential"     : self.use_symlog_for_potential,
            "use_symlog_for_energy"        : self.use_symlog_for_energy,
            "linthresh_for_features"       : self.linthresh_for_features,
            "linthresh_for_potential"      : self.linthresh_for_potential,
            "linthresh_for_energy"         : self.linthresh_for_energy,
            "exclude_rho_from_nn"          : self.exclude_rho_from_nn,
            "loss_norm_for_potential"      : self.loss_norm_for_potential,
            "loss_norm_for_energy"         : self.loss_norm_for_energy,

            # scaler summaries (json friendly)
            "scaler_features_summary"      : self._scaler_summary(self.scaler_for_features),
            "scaler_for_potential_summary" : self._scaler_summary(self.scaler_for_potential),
            "scaler_for_energy_summary"    : self._scaler_summary(self.scaler_for_energy),
        }


    @staticmethod
    def _jsonable(value: Any) -> Any:
        if isinstance(value, (np.integer, np.floating)):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, dict):
            return {k: MLXCCalculator._jsonable(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [MLXCCalculator._jsonable(v) for v in value]
        return value


    @staticmethod
    def _scaler_summary(scaler: Any) -> Optional[Dict[str, Any]]:
        if scaler is None:
            return None
        summary = {
            "type": type(scaler).__name__,
        }
        # common sklearn scaler attributes
        for key in ("mean_", "scale_", "var_", "center_", "quantile_range_"):
            if hasattr(scaler, key):
                summary[key] = MLXCCalculator._jsonable(getattr(scaler, key))
        return summary


    def _resolve_scalers_path(self, model_dir: Union[str, Path]) -> Path:
        model_dir = Path(model_dir)
        return model_dir / f"{self.model.model_name}_scalers.pkl"


    def _resolve_metadata_path(self, model_dir: Union[str, Path]) -> Path:
        model_dir = Path(model_dir)
        return model_dir / f"{self.model.model_name}_metadata.json"


    def _resolve_model_paths(self, model_dir: Union[str, Path]) -> List[Path]:
        model_dir = Path(model_dir)
        model_file = model_dir / f"{self.model.model_name}_model.{self.model.weights_ext}"
        config_file = model_dir / f"{self.model.model_name}_config.{self.model.config_ext}"
        scalers_file = self._resolve_scalers_path(model_dir)
        metadata_file = self._resolve_metadata_path(model_dir)
        return [model_file, config_file, scalers_file, metadata_file]


    @property
    def model_kind(self) -> str:
        return self.model.model_kind 