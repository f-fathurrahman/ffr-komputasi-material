"""
Base class for Delta Learning XC Functional Corrections

This module provides a base class for integrating neural network models
into XC functional evaluations. Users should inherit from DeltaXCEvaluator
and implement the required methods to integrate their trained models.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional, Union, Literal, Type, List, TYPE_CHECKING
import numpy as np

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

try:
    import jax
    JAX_AVAILABLE = True
except ImportError:
    JAX_AVAILABLE = False


# Lazy import to avoid circular: data_manager -> data_generation -> scf -> ml_xc -> base -> data_manager
if TYPE_CHECKING:
    from ...data.data_manager import VxcDataLoader, ExcDataLoader


from ...data.data_loading import (
    DataLoader,
    NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL,
    NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY,
    format_invalid_feature_error,
)

# Error messages for XCBaseModel
MODEL_KIND_NOT_STRING_ERROR = \
    "parameter 'model_kind' must be a string, get {} instead."
MODEL_KIND_NOT_VALID_ERROR = \
    "parameter 'model_kind' must be 'potential' or 'energy', get {} instead."

FEATURES_LIST_NOT_LIST_ERROR = \
    "parameter 'features_list' must be a list, get {} instead."
FEATURE_NOT_STRING_ERROR = \
    "parameter 'features_list' must be a list of strings, get {} instead."

WEIGHTS_EXT_NOT_STRING_ERROR = \
    "parameter 'weights_ext' must be a string, get {} instead."
CONFIG_EXT_NOT_STRING_ERROR = \
    "parameter 'config_ext' must be a string, get {} instead."
MODEL_INSTANCE_AND_CLS_PROVIDED_ERROR = \
    "provide either model_instance or model_cls/model_init_kwargs, not both."
MODEL_NAME_REQUIRED_ERROR = \
    "parameter 'model_name' is required when loading from disk."
MODEL_INIT_KWARGS_NOT_PROVIDED_ERROR = \
    "model_init_kwargs is required when model_cls is provided."
MODEL_INIT_KWARGS_NOT_DICT_ERROR = \
    "model_init_kwargs must be a dictionary, get {} instead."
MODEL_INITIALIZATION_ERROR = \
    "Error initializing model from cls and init_kwargs: {}."
MODEL_NAME_REQUIRED_ERROR = \
    "parameter 'model_name' is required when loading from disk."
MODEL_NAME_NOT_STRING_ERROR = \
    "parameter 'model_name' must be a string, get {} instead."
DEVICE_NOT_STRING_ERROR = \
    "parameter 'device' must be a string, get {} instead."
WEIGHTS_FILE_NOT_FOUND_ERROR = \
    "model's weights file {} does not exist."
CONFIG_FILE_NOT_FOUND_ERROR = \
    "model's config file {} does not exist."

# Training parameters error messages
EPOCHS_NOT_INTEGER_ERROR = \
    "parameter 'epochs' must be an integer, get {} instead."
EPOCHS_NOT_GREATER_THAN_ZERO_INTEGER_ERROR = \
    "parameter 'epochs' must be an integer greater than 0, get {} instead."
LR_NOT_FLOAT_ERROR = \
    "parameter 'lr' must be a float, get {} instead."
LR_NOT_GREATER_THAN_ZERO_FLOAT_ERROR = \
    "parameter 'lr' must be a float greater than 0, get {} instead."
PATIENCE_NOT_INTEGER_ERROR = \
    "parameter 'patience' must be an integer, get {} instead."
PATIENCE_NOT_GREATER_THAN_ZERO_INTEGER_ERROR = \
    "parameter 'patience' must be an integer greater than 0, get {} instead."

BATCH_SIZE_NOT_INTEGER_ERROR = \
    "parameter 'batch_size' must be an integer, get {} instead."
BATCH_SIZE_NOT_POSITIVE_INTEGER_ERROR = \
    "parameter 'batch_size' must be an integer greater than 0, get {} instead."

SHUFFLE_NOT_BOOL_ERROR = \
    "parameter 'shuffle' must be a boolean, get {} instead."





UNKNOWN_INITIALIZATION_METHOD_ERROR = \
    "Unknown initialization method. Use exactly one of: " \
    "(1) construct from class    : model_cls + model_init_kwargs; " \
    "(2) construct from instance : model_instance; " \
    "(3) load from disk          : model_dir + model_name."


ModelKind   = Literal["potential", "energy"]
ConfigExt   = Literal["json", "yaml", "yml"]
WeightsExt  = Literal["pth", "pt", "pkl", "joblib"]
ModelSource = Literal["from_cls", "from_instance", "from_disk"]
DataLoaderType = Union["VxcDataLoader", "ExcDataLoader"]



@dataclass(frozen=True)
class ModelIOFiles:
    """
    Two-file IO layout for a model.
    """
    weights_path: Path
    config_path: Path




@dataclass(frozen=True)
class EvaluationMetrics:
    """
    Evaluation metrics for the model.
    """
    mae  : float
    rmse : float
    r2   : float
    max_relative_error : float


    def print_info(self, model_name: str = "ML model", target_name: str = "") -> None:
        """Print metrics in a readable format.
        
        Parameters
        ----------
        metrics : dict
            Dictionary containing metrics (mae, rmse, r2, max_relative_error)
        model_name : str
            Name of the model
        target_name : str
            Name of the target variable
        """
        info_length = 75
        print(f"{'=' * info_length}")
        print(f"Metrics for {model_name}".center(info_length))
        print(f"{'=' * info_length}")
        
        if target_name != "":
            print(f"Target ({target_name}):".center(info_length))
        
        print(f"  MAE:                {self.mae:.6f}".center(info_length))
        print(f"  RMSE:               {self.rmse:.6f}".center(info_length))
        print(f"  R²:                 {self.r2:.6f}".center(info_length))
        print(f"  Max Relative Error: {self.max_relative_error:.6f}".center(info_length))
        print(f"{'=' * info_length}\n")


@dataclass(frozen=True)
class PotentialEvaluationMetrics:
    """
    Evaluation metrics for potential model.
    """
    mae  : float
    rmse : float
    r2   : float
    max_relative_error : float


    def print_info(self, model_name: str = "ML model", target_name: str = "") -> None:
        """Print metrics in a readable format.
        
        Parameters
        ----------
        model_name : str
            Name of the model
        target_name : str
            Name of the target variable
        """
        info_length = 75
        print(f"{'=' * info_length}")
        print(f"Metrics for {model_name}".center(info_length))
        print(f"{'=' * info_length}")
        
        if target_name != "":
            print(f"Target ({target_name}):".center(info_length))
        
        print(f"  MAE:                {self.mae:.6f}".center(info_length))
        print(f"  RMSE:               {self.rmse:.6f}".center(info_length))
        print(f"  R²:                 {self.r2:.6f}".center(info_length))
        print(f"  Max Relative Error: {self.max_relative_error:.6f}".center(info_length))
        print(f"{'=' * info_length}\n")


@dataclass(frozen=True)
class EnergyEvaluationMetrics:
    """
    Evaluation metrics for energy model, including potential, energy, and total metrics.
    """
    energy   : EvaluationMetrics
    total    : EvaluationMetrics
    potential: Optional[EvaluationMetrics] = None


    def print_info(self, model_name: str = "ML model", include_potential: bool = True) -> None:
        """Print metrics in a table format.
        
        Parameters
        ----------
        model_name : str
            Name of the model
        include_potential : bool
            Whether potential metrics are included
        """
        info_length = 100
        
        # Header
        print(f"{'=' * info_length}")
        print(f"Evaluation Metrics for {model_name}".center(info_length))
        print(f"{'=' * info_length}")
        
        # Table header
        header = f"{'Metric':<25} {'Potential':<25} {'Energy':<25} {'Total':<25}"
        print(header)
        print("-" * info_length)
        
        # Format function for values
        def format_value(metrics: Optional[EvaluationMetrics], attr: str) -> str:
            if metrics is None:
                return "None"
            value = getattr(metrics, attr)
            if value == float('inf'):
                return "Inf"
            return f"{value:.6e}"
        
        # Print each metric row
        print(f"{'MAE':<25} {format_value(self.potential, 'mae'):<25} {format_value(self.energy, 'mae'):<25} {format_value(self.total, 'mae'):<25}")
        print(f"{'RMSE':<25} {format_value(self.potential, 'rmse'):<25} {format_value(self.energy, 'rmse'):<25} {format_value(self.total, 'rmse'):<25}")
        print(f"{'R²':<25} {format_value(self.potential, 'r2'):<25} {format_value(self.energy, 'r2'):<25} {format_value(self.total, 'r2'):<25}")
        print(f"{'Max Relative Error':<25} {format_value(self.potential, 'max_relative_error'):<25} {format_value(self.energy, 'max_relative_error'):<25} {format_value(self.total, 'max_relative_error'):<25}")
        
        print(f"{'=' * info_length}\n")



class XCBaseModel(ABC):
    """
    Base model with three initialization paths:
    1) model_cls + model_init_kwargs -> build an instance
    2) model_instance -> use provided instance directly
    3) model_dir + model_name -> load from disk (weights + config)

    The chosen path is recorded in self.model_source.
    """

    model_kind: ModelKind = "potential"

    def __init__(
        self,
        model_kind        : ModelKind,
        features_list     : List[str],
        weights_ext       : WeightsExt,
        config_ext        : ConfigExt,
        model_cls         : Optional[Type[Any]]      = None,  # construct from cls
        model_init_kwargs : Optional[Dict[str, Any]] = None,  # construct from cls
        model_instance    : Optional[Any]            = None,  # construct from instance
        model_name        : Optional[str]            = None,  # load from disk
        model_dir         : Optional[str]            = None,  # load from disk
        device            : Optional[str]            = None,  # load from disk
    ) -> None:

        # Type and value checks
        self.check_model_kind(model_kind)
        features_list = self.check_and_normalize_features_list(features_list, model_kind)
        self.check_weights_ext(weights_ext)
        self.check_config_ext(config_ext)
        self.check_model_init_params(model_cls, model_init_kwargs, model_instance, model_name, model_dir)
        self.check_model_name(model_name)
        self.check_device(device)

        # First method: construct from class
        if model_cls is not None:
            model_instance = model_cls(**model_init_kwargs)
            self.model_source = "from_cls"
        
        # Second method: construct from instance
        elif model_instance is not None:
            model_cls = model_instance.__class__
            self.model_source = "from_instance"
            raise NotImplementedError("Check this branch later")
        
        # Third method: load from disk
        elif model_dir is not None:
            files = self._resolve_io_files(model_dir, weights_ext=weights_ext, config_ext=config_ext)
            assert files.weights_path.exists(), \
                WEIGHTS_FILE_NOT_FOUND_ERROR.format(files.weights_path)
            assert files.config_path.exists(), \
                CONFIG_FILE_NOT_FOUND_ERROR.format(files.config_path)
            model_instance = self.load_model(model_dir, model_name=model_name)
            self.model_source = "from_disk"
            raise NotImplementedError("Check this branch later")

        # unknown initialization method
        else:
            raise ValueError(UNKNOWN_INITIALIZATION_METHOD_ERROR)

        # set model_name if not provided
        if model_name is None and model_cls is not None and hasattr(model_cls, "__name__"):
            model_name = model_cls.__name__

        # Set attributes
        self.model             = model_instance
        self.model_kind        = model_kind
        self.features_list     = features_list

        self.weights_ext       = weights_ext
        self.config_ext        = config_ext
        self.model_cls         = model_cls
        self.model_init_kwargs = model_init_kwargs
        self.model_name        = model_name if model_name is not None else "xc"
        self.model_dir         = Path(model_dir)
        self.device            = device




    @staticmethod
    def check_model_kind(model_kind) -> None:
        assert isinstance(model_kind, str), \
            MODEL_KIND_NOT_STRING_ERROR.format(type(model_kind))
        if model_kind not in ("potential", "energy"):
            raise ValueError(MODEL_KIND_NOT_VALID_ERROR.format(model_kind))


    @staticmethod
    def check_and_normalize_features_list(features_list, model_kind) -> List[str]:
        assert isinstance(features_list, list), \
            FEATURES_LIST_NOT_LIST_ERROR.format(type(features_list))
        for feature in features_list:
            assert isinstance(feature, str), \
                FEATURE_NOT_STRING_ERROR.format(type(feature))

        normalized_features = DataLoader.check_and_normalize_features_list(features_list)
        
        if model_kind == "potential":
            for feature in normalized_features:
                if feature not in NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL:
                    raise ValueError(
                        format_invalid_feature_error(
                            feature,
                            NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL,
                        )
                    )
        elif model_kind == "energy":
            for feature in normalized_features:
                if feature not in NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY:
                    raise ValueError(
                        format_invalid_feature_error(
                            feature,
                            NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY,
                        )
                    )
        else:
            raise ValueError(MODEL_KIND_NOT_VALID_ERROR.format(model_kind))
        return normalized_features



    @staticmethod
    def check_weights_ext(weights_ext) -> None:
        assert isinstance(weights_ext, str), \
            WEIGHTS_EXT_NOT_STRING_ERROR.format(type(weights_ext))
        
        # TODO: check if the weights extension is valid
        pass


    @staticmethod
    def check_config_ext(config_ext) -> None:
        assert isinstance(config_ext, str), \
            CONFIG_EXT_NOT_STRING_ERROR.format(type(config_ext))
        
        # TODO: check if the config extension is valid
        pass



    @staticmethod
    def check_model_init_params(
        model_cls         : Optional[Type[Any]]      = None,
        model_init_kwargs : Optional[Dict[str, Any]] = None,
        model_instance    : Optional[Any]            = None,
        model_name        : Optional[str]            = None,
        model_dir         : Optional[str]            = None,
    ) -> None:

        # check if the model is initialized in the correct way
        if model_instance is not None and (model_cls is not None or model_init_kwargs is not None):
            raise ValueError(MODEL_INSTANCE_AND_CLS_PROVIDED_ERROR)
        if model_instance is None and model_cls is None and model_name is None:
            raise ValueError(MODEL_NAME_REQUIRED_ERROR)
        
        # First method: from cls and init_kwargs
        if model_cls is not None:
            assert model_init_kwargs is not None, \
                MODEL_INIT_KWARGS_NOT_PROVIDED_ERROR
            assert isinstance(model_init_kwargs, dict), \
                MODEL_INIT_KWARGS_NOT_DICT_ERROR.format(type(model_init_kwargs))
            try:
                model_cls(**model_init_kwargs)
            except Exception as e:
                raise ValueError(MODEL_INITIALIZATION_ERROR.format(e))
        
        # Second method: from instance
        elif model_instance is not None:
            assert model_cls is None and model_init_kwargs is None, \
                MODEL_INSTANCE_AND_CLS_PROVIDED_ERROR

        # Third method: from disk
        elif model_dir is not None:
            assert model_name is not None, \
                MODEL_NAME_REQUIRED_ERROR
        
        # unknown initialization method
        else:
            raise ValueError(UNKNOWN_INITIALIZATION_METHOD_ERROR)


    @staticmethod
    def check_model_name(model_name) -> None:
        if model_name is not None:
            assert isinstance(model_name, str), \
                MODEL_NAME_NOT_STRING_ERROR.format(type(model_name))


    @staticmethod
    def check_device(device) -> None:
        assert isinstance(device, str), \
            DEVICE_NOT_STRING_ERROR.format(type(device))
        
        # TODO: check if the device is valid
        pass


    @staticmethod
    def check_epochs(epochs) -> None:
        assert isinstance(epochs, int), \
            EPOCHS_NOT_INTEGER_ERROR.format(type(epochs))
        assert epochs > 0, \
            EPOCHS_NOT_GREATER_THAN_ZERO_INTEGER_ERROR.format(epochs)


    @staticmethod
    def check_lr(lr) -> None:
        assert isinstance(lr, float), \
            LR_NOT_FLOAT_ERROR.format(type(lr))
        assert lr > 0, \
            LR_NOT_GREATER_THAN_ZERO_FLOAT_ERROR.format(lr)


    @staticmethod
    def check_patience(patience) -> None:
        assert isinstance(patience, int), \
            PATIENCE_NOT_INTEGER_ERROR.format(type(patience))
        assert patience > 0, \
            PATIENCE_NOT_GREATER_THAN_ZERO_INTEGER_ERROR.format(patience)


    @staticmethod
    def check_batch_size(batch_size) -> None:
        assert isinstance(batch_size, int), \
            BATCH_SIZE_NOT_INTEGER_ERROR.format(type(batch_size))
        assert batch_size > 0, \
            BATCH_SIZE_NOT_POSITIVE_INTEGER_ERROR.format(batch_size)


    @staticmethod
    def check_shuffle(shuffle) -> None:
        assert isinstance(shuffle, bool), \
            SHUFFLE_NOT_BOOL_ERROR.format(type(shuffle))


    def _resolve_io_files(
        self,
        model_dir   : Optional[Union[str, Path]] = None,
        weights_ext : Optional[WeightsExt]       = None,
        config_ext  : Optional[ConfigExt]        = None,
    ) -> ModelIOFiles:
        """
        Resolve the IO files for the model.
        """
        model_dir   = Path(model_dir) if model_dir   is not None else self.model_dir
        weights_ext = weights_ext     if weights_ext is not None else self.weights_ext
        config_ext  = config_ext      if config_ext  is not None else self.config_ext
        weights_path = model_dir / f"{self.model_name}.{weights_ext}"
        config_path  = model_dir / f"{self.model_name}.{config_ext}"
        return ModelIOFiles(weights_path=weights_path, config_path=config_path)


    def model_exists(
        self,
        model_dir   : Union[str, Path],
        weights_ext : Optional[WeightsExt] = None,
        config_ext  : Optional[ConfigExt]  = None,
    ) -> bool:
        """
        Check if the model files exist.
        """
        files = self._resolve_io_files(model_dir, weights_ext=weights_ext, config_ext=config_ext)
        return files.weights_path.exists() and files.config_path.exists()



    @abstractmethod
    def save_model(
        self,
        model_dir: Union[str, Path],
        *,
        config: Dict[str, Any],
        config_ext: Optional[ConfigExt] = None,
    ) -> ModelIOFiles:
        """
        Save model weights and a single config file (json/yaml).
        """
        raise NotImplementedError


    # @abstractmethod
    def load_model(
        self,
        model_dir  : Union[str, Path],
        model_name : Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Load model weights and config from the two-file layout.
        """
        raise NotImplementedError



    # @abstractmethod
    def eval_model(
        self,
        data_loader: Union["VxcDataLoader", "ExcDataLoader"],
    ) -> EvaluationMetrics:
        """
        Evaluate model on a data loader and return metrics.
        """
        raise NotImplementedError


    @abstractmethod
    def train(
        self,
        train_loader : DataLoaderType,
        val_loader   : Optional[DataLoaderType] = None,
        epochs       : int   = 100,
        lr           : float = 1e-3,
        patience     : int   = 100,
    ) -> Dict[str, Any]:
        """
        Train the model and return a training summary (loss curves, best epoch, etc.).
        """
        self.check_epochs(epochs)
        self.check_lr(lr)
        self.check_patience(patience)

        raise NotImplementedError


    @abstractmethod
    def forward(self, features: np.ndarray) -> np.ndarray:
        """
        Forward pass of the model.
        """
        raise NotImplementedError

