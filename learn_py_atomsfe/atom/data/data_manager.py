
__author__ = "Qihao Cheng"

import os
import traceback
import warnings
import numpy as np
from typing import Optional, Tuple, List, Dict, Any, Union, Literal, get_args, get_origin
from dataclasses import dataclass, field

# Import from sklearn
try:
    from sklearn.preprocessing import StandardScaler, RobustScaler
    from sklearn.model_selection import train_test_split
    SKLEARN_AVAILABLE = True
    ScalerType = Optional[Union[StandardScaler, RobustScaler]]

except ImportError:
    SKLEARN_AVAILABLE = False
    ScalerType = Any

# Error messages for AtomicDataset
DATA_ROOT_NOT_STRING_ERROR = \
    "parameter 'data_root' must be a string, get type {} instead."
DATA_ROOT_NOT_EXIST_ERROR = \
    "data root '{}' does not exist. Please create an empty directory if you want to generate data, or check if the directory is correct."
SCF_XC_FUNCTIONAL_NOT_STRING_ERROR = \
    "parameter 'scf_xc_functional' must be a string, get type {} instead."
SCF_XC_FUNCTIONAL_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'scf_xc_functional' must be in {}, get {} instead."
FORWARD_PASS_XC_FUNCTIONAL_NOT_NONE_STR_OR_LIST_ERROR = \
    "parameter 'forward_pass_xc_functionals' must be None, string, or list of strings, get type {} instead."
FORWARD_PASS_XC_FUNCTIONAL_LIST_NOT_LIST_ERROR = \
    "parameter 'forward_pass_xc_functional_list' must be a list, get {} instead."
FORWARD_PASS_XC_FUNCTIONAL_LIST_NOT_LIST_OF_STRINGS_ERROR = \
    "parameter 'forward_pass_xc_functional_list' must be a list of strings, get {} instead."
FORWARD_PASS_XC_FUNCTIONAL_LIST_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'forward_pass_xc_functional_list' must be in {}, get {} instead."
FEATURES_LIST_NOT_LIST_ERROR = \
    "parameter 'features_list' must be a list, get {} instead."
FEATURES_LIST_NOT_LIST_OF_STRINGS_ERROR = \
    "parameter 'features_list' must be a list of strings, get {} instead."
FEATURES_LIST_MUST_CONTAIN_RHO_ERROR = \
    "parameter 'features_list' must contain 'rho' when initializing ComputationConfig."


RADIUS_CUTOFF_RHO_THRESHOLD_NOT_FLOAT_ERROR = \
    "parameter 'radius_cutoff_rho_threshold' must be a float, get type {} instead."
RADIUS_CUTOFF_V_X_THRESHOLD_NOT_FLOAT_ERROR = \
    "parameter 'radius_cutoff_v_x_threshold' must be a float, get type {} instead."
RADIUS_CUTOFF_V_C_THRESHOLD_NOT_FLOAT_ERROR = \
    "parameter 'radius_cutoff_v_c_threshold' must be a float, get type {} instead."
RADIUS_CUTOFF_RHO_THRESHOLD_NOT_POSITIVE_ERROR = \
    "parameter 'radius_cutoff_rho_threshold' must be positive, get {} instead."
RADIUS_CUTOFF_V_X_THRESHOLD_NOT_POSITIVE_ERROR = \
    "parameter 'radius_cutoff_v_x_threshold' must be positive, get {} instead."
RADIUS_CUTOFF_V_C_THRESHOLD_NOT_POSITIVE_ERROR = \
    "parameter 'radius_cutoff_v_c_threshold' must be positive, get {} instead."

SMOOTH_RADIUS_THRESHOLD_NOT_FLOAT_ERROR = \
    "parameter 'smooth_radius_threshold' must be a float, get type {} instead."
SMOOTH_RADIUS_THRESHOLD_NOT_POSITIVE_ERROR = \
    "parameter 'smooth_radius_threshold' must be positive, get {} instead."
SMOOTH_METHOD_NOT_STRING_ERROR = \
    "parameter 'smooth_method' must be a string, get type {} instead."
SMOOTH_METHOD_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'smooth_method' must be in {}, get {} instead."
SMOOTH_KWARGS_NOT_DICT_ERROR = \
    "parameter 'smooth_kwargs' must be a dict, get type {} instead."


CONFIGURATION_DATA_LIST_NOT_LIST_ERROR = \
    "parameter 'configuration_data_list' must be a list, get {} instead."
CONFIGURATION_DATA_LIST_NOT_LIST_OF_SINGLE_CONFIGURATION_DATA_ERROR = \
    "parameter 'configuration_data_list' must be a list of SingleConfigurationData, invalid element types: {}."

ATOMIC_NUMBER_NOT_INTEGER_ERROR = \
    "atomic number must be an integer, get {} instead."
ATOMIC_NUMBER_NOT_IN_DATASET_ERROR = \
    "atomic number {} is not in the dataset."
FUNCTIONAL_NOT_IN_DATASET_ERROR = \
    "functional '{}' is not in dataset: scf={}, forward_pass={}."

FEATURE_NAME_NOT_STRING_ERROR = \
    "parameter 'feature_name' must be a string, get {} instead."
FEATURE_NAME_NOT_IN_FEATURES_LIST_ERROR = \
    "parameter 'feature_name' must be in {}, get {} instead."

WEIGHTS_MODE_NOT_STRING_ERROR = \
    "parameter 'weights_mode' must be a string, get type {} instead."
WEIGHTS_MODE_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'weights_mode' must be in {}, get {} instead."

TARGET_COMPONENT_NOT_PROVIDED_ERROR = \
    "parameter 'target_component' must be provided for symlog weights mode, get None instead."
TARGET_FUNCTIONAL_NOT_PROVIDED_ERROR = \
    "parameter 'target_functional' must be provided for symlog weights mode, get None instead."
LINTHRESH_FOR_TARGETS_NOT_PROVIDED_ERROR = \
    "parameter 'linthresh_for_targets' must be provided for symlog weights mode, get None instead."
LOSS_NORM_NOT_IN_VALID_LIST_ERROR = \
    "parameter '{}' must be in {}, get {} instead."

# Error messages for DataManager
ATOMIC_NUMBER_LIST_NOT_LIST_ERROR = \
    "parameter 'atomic_number_list' must be a list, get {} instead."
ATOMIC_NUMBER_LIST_NOT_LIST_OF_INTEGERS_ERROR = \
    "parameter 'atomic_number_list' must be a list of integers, get {} instead."
ATOMIC_NUMBER_LIST_NOT_IN_VALID_RANGE_ERROR = \
    "parameter 'atomic_number_list' must be a list of integers in the range 1-92, get {} instead."
ATOMIC_NUMBER_LIST_AND_CONFIGURATION_INDEX_LIST_BOTH_SPECIFIED_ERROR = \
    "Cannot specify both 'atomic_number_list' (deprecated) and 'configuration_index_list'. " \
    "Use 'configuration_index_list' only."
USE_RADIUS_CUTOFF_NOT_BOOL_ERROR = \
    "parameter 'use_radius_cutoff' must be a boolean, get {} instead."
USE_FEATURE_ROUND_OFF_NOT_BOOL_ERROR = \
    "parameter 'use_feature_round_off' must be a boolean, get {} instead."
SMOOTH_XC_DATA_NOT_BOOL_ERROR = \
    "parameter 'smooth_xc_data' must be a boolean, get {} instead."
CLOSE_SHELL_ONLY_NOT_BOOL_ERROR = \
    "parameter 'close_shell_only' must be a boolean, get {} instead."
SPIN_UNPOLARIZED_ONLY_NOT_BOOL_ERROR = \
    "parameter 'spin_unpolarized_only' must be a boolean, get {} instead."
INCLUDE_ENERGY_DENSITY_NOT_BOOL_ERROR = \
    "parameter 'include_energy_density' must be a boolean, get {} instead."
INCLUDE_INTERMEDIATE_NOT_BOOL_ERROR = \
    "parameter 'include_intermediate' must be a boolean, get {} instead."
SMOOTH_REDUCED_LAPLACIAN_NOT_BOOL_ERROR = \
    "parameter 'smooth_reduced_laplacian' must be a boolean, get {} instead."
SMOOTH_RHO_FOR_DERIVATIVES_NOT_BOOL_ERROR = \
    "parameter 'smooth_rho_for_derivatives' must be a boolean, get {} instead."
PRINT_DEBUG_INFO_NOT_BOOL_ERROR = \
    "parameter 'print_debug_info' must be a boolean, get {} instead."
PRINT_SUMMARY_NOT_BOOL_ERROR = \
    "parameter 'print_summary' must be a boolean, get {} instead."

REFERENCE_FUNCTIONAL_NOT_STRING_ERROR = \
    "parameter 'reference_functional' must be a string, get {} instead."
REFERENCE_FUNCTIONAL_NOT_IN_DATASET_ERROR = \
    "parameter 'reference_functional' '{}' is not in dataset: scf={}, forward_pass={}."
TARGET_FUNCTIONAL_NOT_STRING_ERROR = \
    "parameter 'target_functional' must be a string, get {} instead."
TARGET_FUNCTIONAL_NOT_IN_DATASET_ERROR = \
    "parameter 'target_functional' '{}' is not in dataset: scf={}, forward_pass={}."
TARGET_COMPONENT_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'target_component' must be in {}, get {} instead."
PREPARE_ENERGY_DATALOADER_REQUIRES_INCLUDE_ENERGY_DENSITY_ERROR = \
    "prepare_energy_dataloader requires include_energy_density=True when loading data."
SPLIT_EMPTY_VAL_OR_TEST_ERROR = \
    "Split produced empty val (n={}) or test (n={}) set. Use fewer atoms in ensure_train_atoms or add more data."
SKLEARN_NOT_AVAILABLE_FOR_DATA_PREPROCESSING_ERROR = \
    "sklearn is required for data preprocessing, but it is not available. To install sklearn, run 'pip install scikit-learn'."


# Warning messages
V_XC_IS_ALREADY_SMOOTHED_WARNING = \
    "WARNING: V_xc is already smoothed for this dataset, skipping smoothing"
E_XC_IS_ALREADY_SMOOTHED_WARNING = \
    "WARNING: E_xc is already smoothed for this dataset, skipping smoothing"
POTENTIAL_WEIGHTS_DATA_ALREADY_UPDATED_WARNING = \
    "WARNING: Potential weights data is already updated for this dataset, update again will overwrite the existing weights data."
ENERGY_WEIGHTS_DATA_ALREADY_UPDATED_WARNING = \
    "WARNING: Energy weights data is already updated for this dataset, update again will overwrite the existing weights data."
SMOOTH_REDUCED_LAPLACIAN_NO_LAP_RHO_REDUCED_WARNING = \
    "WARNING: smooth_reduced_laplacian=True has no effect when 'lap_rho_reduced' is not in features_list."
SMOOTH_RHO_FOR_DERIVATIVES_NO_DERIV_FEATURES_WARNING = \
    "WARNING: smooth_rho_for_derivatives=True only changes rho and derivative features; " \
    "no grad/lap features in features_list (rho is still asymptotically smoothed)."


# Deprecated arguments warning messages
USE_CUTOFF_DEPRECATED_WARNING = \
    "WARNING: parameter 'use_cutoff' is now deprecated, use parameter 'use_radius_cutoff' instead"
CUTOFF_RHO_THRESHOLD_DEPRECATED_WARNING = \
    "WARNING: parameter 'cutoff_rho_threshold' is now deprecated, use parameter 'radius_cutoff_rho_threshold' instead"
CUTOFF_V_X_THRESHOLD_DEPRECATED_WARNING = \
    "WARNING: parameter 'cutoff_v_x_threshold' is now deprecated, use parameter 'radius_cutoff_v_x_threshold' instead"
CUTOFF_V_C_THRESHOLD_DEPRECATED_WARNING = \
    "WARNING: parameter 'cutoff_v_c_threshold' is now deprecated, use parameter 'radius_cutoff_v_c_threshold' instead"
FINITE_ELEMENTS_DEPRECATED_WARNING = \
    "WARNING: parameter 'finite_elements' is now deprecated, use 'finite_elements_number' instead."
FINITE_ELEMENTS_NUMBER_AND_FINITE_ELEMENTS_BOTH_SPECIFIED_ERROR = \
    "Cannot specify both 'finite_elements_number' and deprecated 'finite_elements'. Use 'finite_elements_number' only."
ATOMIC_NUMBER_LIST_DEPRECATED_WARNING = \
    "WARNING: parameter 'atomic_number_list' is deprecated. Use 'configuration_index_list' instead. \n" \
    "\t Atomic numbers are now read from meta.json files, not from folder names."

SCALE_TARGETS_DEPRECATED_WARNING = \
    "WARNING: parameter 'scale_targets' is now deprecated, use 'scale_potential' instead."
SCALER_TYPE_TARGETS_DEPRECATED_WARNING = \
    "WARNING: parameter 'scaler_type_targets' is now deprecated, use 'scaler_type_for_potential' instead."
USE_SYMLOG_TARGETS_DEPRECATED_WARNING = \
    "WARNING: parameter 'use_symlog_targets' is now deprecated, use 'use_symlog_for_potential' instead."
SCALER_TYPE_FEATURES_DEPRECATED_WARNING = \
    "WARNING: parameter 'scaler_type_features' is now deprecated, use 'scaler_type_for_features' instead."
SCALER_KWARGS_TARGETS_DEPRECATED_WARNING = \
    "WARNING: parameter 'scaler_kwargs_targets' is now deprecated, use 'scaler_kwargs_for_potential' instead."
USE_SYMLOG_FEATURES_DEPRECATED_WARNING = \
    "WARNING: parameter 'use_symlog_features' is now deprecated, use 'use_symlog_for_features' instead."
LINTHRESH_FEATURES_DEPRECATED_WARNING = \
    "WARNING: parameter 'linthresh_features' is now deprecated, use 'linthresh_for_features' instead."
LINTHRESH_TARGETS_DEPRECATED_WARNING = \
    "WARNING: parameter 'linthresh_targets' is now deprecated, use 'linthresh_for_potential' instead."
MIN_WEIGHT_RATIO_DEPRECATED_WARNING = \
    "WARNING: parameter 'min_weight_ratio' is now deprecated, use 'min_weight_ratio_for_potential' instead."
SMOOTH_VXC_DEPRECATED_WARNING = \
    "WARNING: parameter 'smooth_vxc' is now deprecated, use 'smooth_xc_data' instead."

SPLIT_DATA_BY_ATOM_DATA_LEAKAGE_WARNING = \
    """
    [ATOM WARNING] split_data_by_atom() may cause data leakage.

    This method splits an already-built VxcDataLoader whose scaling/symlog
    was applied to the full dataset. As a result:

        • Scalers (e.g. StandardScaler/RobustScaler) were fit on data that
          includes validation and test samples, so val/test statistics leak
          into preprocessing.
        • Potential weights (and any normalization using them) were also
          computed using the full data.

    Recommended: use prepare_potential_dataloader(..., split_train_val_test=True)
    instead. That path:

        1) Splits by atom first to obtain train/val/test masks,
        2) Fits scalers and weights only on the training set,
        3) Transforms val/test with the same fitted objects (no refit),

    so no information from val/test is used during preprocessing.
    """





def format_error_message(e: Exception, context: str = "") -> Tuple[str, str]:
    """Format error message with type and enhanced context.
    
    Args:
        e: Exception instance
        context: Additional context string
        
    Returns:
        Tuple of (error_summary, full_traceback)
    """
    error_type = type(e).__name__
    error_msg = str(e)
    
    # Enhance KeyError messages
    if isinstance(e, KeyError):
        error_msg = f"Missing key '{error_msg}' (possibly a string formatting issue)"
    
    summary = f"[{error_type}] {error_msg}"
    if context:
        summary = f"{context}: {summary}"
    
    return summary, traceback.format_exc()


# Import other modules after format_error_message is defined
# (data_loading.py needs format_error_message, so it must be defined first)
from .data_generation import DataGenerator, XC_FUNCTIONAL_OEP_DEFAULT
from .data_loading import (
    DataLoader,
    SingleConfigurationData,
    FEATURE_ALIASES,
    NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL,
    VALID_POTENTIAL_TARGET_COMPONENTS,
    VALID_ENERGY_TARGET_COMPONENTS,
    PotentialTargetComponent,
    EnergyTargetComponent,
    format_invalid_feature_error,
)
from .data_processing import DataProcessor, VALID_SMOOTH_METHODS
from ..utils.occupation_states import OccupationInfo

# Type aliases
# (v_x, v_c, e_x, e_c), where e_x and e_c are only available if include_energy_density is True
XCDataType = Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]


VALID_WEIGHTS_MODES = {"symlog", "density"}

# Per-grid sample weights when targets are symlog + scaled (potential and energy; see update_*_weights_data).
# L1norm: weighted L1 / MAE alignment; L2norm: weighted MSE alignment (extra σ/|Symlog'| factor).
LossNorm = Literal["L1norm", "L2norm"]
VALID_LOSS_NORMS = frozenset(get_args(LossNorm))


def _validate_loss_norm(param_name: str, value: Any) -> None:
    if value not in VALID_LOSS_NORMS:
        raise ValueError(LOSS_NORM_NOT_IN_VALID_LIST_ERROR.format(param_name, VALID_LOSS_NORMS, value))


class TypeCheckMixin:
    
    @staticmethod
    def check_type(value, expected_type, name):
        def _raise_type_error():
            raise TypeError(
                f"Attribute {name} must be of type {expected_type}, get {type(value)} instead"
            )

        if expected_type is Any:
            return

        origin = get_origin(expected_type)
        if origin is None:
            if not isinstance(value, expected_type):
                _raise_type_error()
            return

        if origin is Union:
            valid_types = tuple(
                t for t in get_args(expected_type) if t is not type(None)
            )
            if value is None:
                return
            if valid_types and isinstance(value, valid_types):
                return
            _raise_type_error()

        if not isinstance(value, origin):
            _raise_type_error()

    @staticmethod
    def check_dim(value, expected_dim, name):
        if value.ndim != expected_dim:
            raise ValueError(f"Attribute '{name}' must have {expected_dim} dimensions, get {value.ndim} instead")
    
    @staticmethod
    def check_shape(value, expected_shape, name):
        if value.shape != expected_shape:
            raise ValueError(f"Attribute '{name}' must have shape {expected_shape}, get {value.shape} instead")
    

    @staticmethod
    def check_is_not_none(value, name, condition_name, condition_value):
        if value is None:
            raise ValueError(f"Attribute '{name}' must NOT be None when '{condition_name}' is {condition_value}, get {value} instead")


@dataclass
class VxcDataLoader(TypeCheckMixin):
    """
    Data class for V_xc data loader. Here, V_xc is the exchange-correlation potential.
    """

    # data attributes (transformed)
    features                    : np.ndarray  # (n_samples, n_features)
    potential_transformed       : np.ndarray  # (n_samples, n_targets)
    weights_for_potential       : np.ndarray  # (n_samples, n_targets)
    atomic_numbers_per_sample   : np.ndarray  # (n_samples,)
    n_electrons_per_sample      : np.ndarray  # (n_samples,)
    
    # parameters for documentation
    target_functional           : str
    target_component            : PotentialTargetComponent
    target_mode                 : str
    reference_functional        : Optional[str]
    features_list               : List[str] = field(default_factory=list)

    # optional parameters for scaling
    scale_features              : bool = True
    scale_potential             : bool = True
    scaler_type_for_features    : str = 'robust'
    scaler_type_for_potential   : str = 'robust'
    scaler_kwargs_for_features  : Dict[str, Any] = field(default_factory=dict)
    scaler_kwargs_for_potential : Dict[str, Any] = field(default_factory=dict)
    scaler_for_features         : ScalerType = None
    scaler_for_potential        : ScalerType = None

    # optional parameters for symlog transformation
    use_symlog_for_features     : bool = True
    use_symlog_for_potential    : bool = True
    linthresh_for_features      : Optional[float] = 0.002
    linthresh_for_potential     : Optional[float] = 0.002

    # symlog+scale sample weights: L1norm (MAE) vs L2norm (MSE) quadrature matching
    loss_norm_for_potential : LossNorm = "L1norm"


    def __post_init__(self):
        # type checks
        self.check_type(self.features                   , np.ndarray      , "features")
        self.check_type(self.potential_transformed      , np.ndarray      , "potential_transformed")
        self.check_type(self.weights_for_potential      , np.ndarray      , "weights_for_potential")
        self.check_type(self.atomic_numbers_per_sample  , np.ndarray      , "atomic_numbers_per_sample")
        self.check_type(self.n_electrons_per_sample     , np.ndarray      , "n_electrons_per_sample")
        self.check_type(self.target_functional          , str             , "target_functional")
        self.check_type(self.target_component           , str             , "target_component")
        self.check_type(self.target_mode                , str             , "target_mode")
        self.check_type(self.reference_functional       , Optional[str]   , "reference_functional")
        self.check_type(self.features_list              , List[str]       , "features_list")
        self.check_type(self.scale_features             , bool            , "scale_features")
        self.check_type(self.scale_potential            , bool            , "scale_potential")
        self.check_type(self.scaler_type_for_features   , str             , "scaler_type_for_features")
        self.check_type(self.scaler_type_for_potential  , str             , "scaler_type_for_potential")
        self.check_type(self.scaler_kwargs_for_features , Dict[str, Any]  , "scaler_kwargs_for_features")
        self.check_type(self.scaler_kwargs_for_potential, Dict[str, Any]  , "scaler_kwargs_for_potential")
        self.check_type(self.scaler_for_features        , ScalerType      , "scaler_for_features")
        self.check_type(self.scaler_for_potential       , ScalerType      , "scaler_for_potential")
        self.check_type(self.use_symlog_for_features    , bool            , "use_symlog_for_features")
        self.check_type(self.use_symlog_for_potential   , bool            , "use_symlog_for_potential")
        self.check_type(self.linthresh_for_features     , Optional[float] , "linthresh_for_features")
        self.check_type(self.linthresh_for_potential    , Optional[float] , "linthresh_for_potential")
        self.check_type(self.loss_norm_for_potential    , str             , "loss_norm_for_potential")

        # dimension checks
        self.check_dim(self.features                  , 2, "features")
        self.check_dim(self.potential_transformed     , 2, "potential_transformed")
        self.check_dim(self.weights_for_potential     , 2, "weights_for_potential")
        self.check_dim(self.atomic_numbers_per_sample , 1, "atomic_numbers_per_sample")
        self.check_dim(self.n_electrons_per_sample    , 1, "n_electrons_per_sample")

        # shape checks
        self.check_shape(self.features             , (self.n_samples, self.n_features), "features")
        self.check_shape(self.potential_transformed, (self.n_samples, self.n_targets) , "potential_transformed")
        self.check_shape(self.weights_for_potential, (self.n_samples, self.n_targets) , "weights_for_potential")
        self.check_shape(self.atomic_numbers_per_sample, (self.n_samples,), "atomic_numbers_per_sample")
        self.check_shape(self.n_electrons_per_sample   , (self.n_samples,), "n_electrons_per_sample")

        # value checks
        if self.scale_features:
            self.check_is_not_none(self.scaler_for_features, "scaler_for_features", "scale_features", True)
        if self.scale_potential:
            self.check_is_not_none(self.scaler_for_potential, "scaler_for_potential", "scale_potential", True)
        if self.use_symlog_for_features:
            self.check_is_not_none(self.linthresh_for_features, "linthresh_for_features", "use_symlog_for_features", True)
        if self.use_symlog_for_potential:
            self.check_is_not_none(self.linthresh_for_potential, "linthresh_for_potential", "use_symlog_for_potential", True)
        _validate_loss_norm("loss_norm_for_potential", self.loss_norm_for_potential)

        # other checks
        assert self.target_mode in ["absolute", "delta"]
        

    def print_info(self, label: str = "Vxc Data Loader"):
        """
        Print information about the V_xc data loader.
        """
        print(f"\n{'='*75}")
        print(f"{label} Summary".center(75))
        print(f"{'='*75}")
        print(f"Number of atomic numbers   : {self.n_atoms}")
        print(f"Number of samples          : {self.n_samples}")
        print(f"Number of features         : {self.n_features}")
        print(f"Number of targets          : {self.n_targets}")
        print(f"Target functional          : {self.target_functional}")
        print(f"Target component           : {self.target_component}")
        print(f"Target mode                : {self.target_mode}")
        print(f"Reference functional       : {self.reference_functional}")
        print(f"Features in features list  : {len(self.features_list)} channels")
        for idx, feature in enumerate(self.features_list):
            suffix = " (repeated)" if feature in self.features_list[:idx] else ""
            print(f"    - Channel {idx + 1}: {feature}{suffix}")
        print()
        print(f"shape of features          : Array of shape {self.features.shape}")
        print(f"shape of potential         : Array of shape {self.potential_transformed.shape}")
        print(f"shape of weights           : Array of shape {self.weights_for_potential.shape}")
        print(f"loss_norm_for_potential    : {self.loss_norm_for_potential}")
        print(f"shape of atomic_numbers    : Array of shape {self.atomic_numbers_per_sample.shape}")
        print(f"shape of n_electrons       : Array of shape {self.n_electrons_per_sample.shape}")
        print(f"{'='*75}")


    def get_features_data(self, feature_name_list: List[str]) -> np.ndarray:
        """
        Get features data for a given feature list.
        """
        for feature_name in feature_name_list:
            assert feature_name in self.features_list, \
                FEATURE_NAME_NOT_IN_FEATURES_LIST_ERROR.format(feature_name, self.features_list)
        
        feature_index_list = [self.features_list.index(feature_name) for feature_name in feature_name_list]
        return self.features[:, feature_index_list]


    @property
    def n_atoms(self) -> int:
        return len(np.unique(self.atomic_numbers_per_sample))

    @property
    def n_features(self) -> int:
        return len(self.features_list)

    @property
    def n_samples(self) -> int:
        return len(self.atomic_numbers_per_sample)

    @property
    def n_targets(self) -> int:
        return self.potential_transformed.shape[1]

    @property
    def atomic_number_list(self) -> List[int]:
        return np.unique(self.atomic_numbers_per_sample).tolist()


    def split_data_by_atom(
        self,
        test_size           : float = 0.2,
        val_size            : float = 0.1,
        random_state        : int   = 42, 
        ensure_train_atoms  : Optional[List[int]] = None
    ) -> Tuple["VxcDataLoader", "VxcDataLoader", "VxcDataLoader"]:
        """
        Split data ensuring atoms don't leak between train/val/test sets.
        
        Parameters
        ----------
        test_size : float
            Proportion of atoms for test set
        val_size : float
            Proportion of atoms for validation set
        random_state : int
            Random seed
        ensure_train_atoms : list or None
            List of atomic numbers that must be in training set (e.g., [0, 1, 2, ..., 20]).

        Returns
        -------
        train_loader, val_loader, test_loader : VxcDataLoader
            Data loaders split by atom without leakage.
        """
        # This method is deprecated and will be removed in the future.
        warnings.warn(SPLIT_DATA_BY_ATOM_DATA_LEAKAGE_WARNING, UserWarning, stacklevel=2)

        atomic_numbers = self.atomic_numbers_per_sample
        n_electrons = self.n_electrons_per_sample
        unique_atoms = np.unique(atomic_numbers)
        
        # Ensure specified atoms are in training set
        if ensure_train_atoms is not None:
            ensure_train_atoms = np.array(ensure_train_atoms)
            # Find atoms that exist in the data and should be in training set
            atoms_guaranteed_train = np.intersect1d(unique_atoms, ensure_train_atoms)
            # Remaining atoms to split
            atoms_to_split = np.setdiff1d(unique_atoms, atoms_guaranteed_train)
            
            if len(atoms_guaranteed_train) > 0:
                print(f"Ensuring atoms {atoms_guaranteed_train.tolist()} are in training set")
        else:
            atoms_guaranteed_train = np.array([], dtype=int)
            atoms_to_split = unique_atoms
        
        # Split remaining atoms into train/val/test
        if len(atoms_to_split) > 0:
            atoms_train_temp, atoms_temp = train_test_split(
                atoms_to_split, test_size=(test_size + val_size), random_state=random_state
            )
            val_ratio = val_size / (test_size + val_size)
            atoms_val, atoms_test = train_test_split(
                atoms_temp, test_size=(1 - val_ratio), random_state=random_state
            )
            
            # Combine guaranteed training atoms with randomly split training atoms
            atoms_train = np.concatenate([atoms_guaranteed_train, atoms_train_temp])
        else:
            # All atoms are guaranteed to be in training set
            atoms_train = atoms_guaranteed_train
            atoms_val = np.array([], dtype=int)
            atoms_test = np.array([], dtype=int)
        
        # Create masks
        train_mask = np.isin(atomic_numbers, atoms_train)
        val_mask   = np.isin(atomic_numbers, atoms_val)
        test_mask  = np.isin(atomic_numbers, atoms_test)

        def _subset(mask: np.ndarray) -> "VxcDataLoader":
            return VxcDataLoader(
                features                    = self.features[mask],
                potential_transformed       = self.potential_transformed[mask],
                weights_for_potential       = self.weights_for_potential[mask],
                atomic_numbers_per_sample   = atomic_numbers[mask],
                n_electrons_per_sample      = n_electrons[mask],
                target_functional           = self.target_functional,
                target_component            = self.target_component,
                target_mode                 = self.target_mode,
                reference_functional        = self.reference_functional,
                features_list               = self.features_list,
                scale_features              = self.scale_features,
                scale_potential             = self.scale_potential,
                scaler_type_for_features    = self.scaler_type_for_features,
                scaler_type_for_potential   = self.scaler_type_for_potential,
                scaler_kwargs_for_features  = self.scaler_kwargs_for_features,
                scaler_kwargs_for_potential = self.scaler_kwargs_for_potential,
                scaler_for_features         = self.scaler_for_features,
                scaler_for_potential        = self.scaler_for_potential,
                use_symlog_for_features     = self.use_symlog_for_features,
                use_symlog_for_potential    = self.use_symlog_for_potential,
                linthresh_for_features      = self.linthresh_for_features,
                linthresh_for_potential     = self.linthresh_for_potential,
                loss_norm_for_potential     = self.loss_norm_for_potential,
            )

        return _subset(train_mask), _subset(val_mask), _subset(test_mask)


@dataclass
class ExcConfiguration(TypeCheckMixin):
    """
    Data class for a single configuration in ExcDataLoader.
    Contains the required fields for energy and optionally potential data.
    
    Note on derivative_matrix:
        If derivative_matrix is None, the shared derivative matrix stored in ExcDataLoader
        (accessible via the parent loader's shared_derivative_matrix attribute) should be used.
        This allows multiple configurations to share the same derivative matrix when they use
        the same grid/basis parameters.
    """
    index                        : int
    configuration_id             : int
    finite_element_number        : int  # n_elem; finite_element_number * quadrature_point_number == n_grid
    quadrature_point_number      : int  # n_quad per element
    quadrature_nodes             : np.ndarray  # (n_grid,)
    quadrature_weights           : np.ndarray  # (n_grid,)
    features                     : np.ndarray  # (n_grid_filtered, n_features) - transformed (symlog + scale)
    features_physical            : np.ndarray  # (n_grid_filtered, n_features) - original untransformed features
    energy_density_transformed   : np.ndarray  # (n_grid_filtered, n_targets)
    weights_for_energy_density   : np.ndarray  # (n_grid_filtered, n_targets)
    total_energy_true            : Optional[np.ndarray] = None  # (n_targets,) integrated E = ∫4πr²w·e; precomputed at load for total_energy loss
    potential_transformed        : Optional[np.ndarray] = None  # (n_grid_filtered, n_targets), optional
    weights_for_potential        : Optional[np.ndarray] = None  # (n_grid_filtered, n_targets), optional
    derivative_matrix            : Optional[np.ndarray] = None  # (n_elem, n_quad, n_quad), optional. If None, use shared_derivative_matrix from ExcDataLoader.


    def __post_init__(self):
        # Type checks
        self.check_type(self.index                      , int                  , "index")
        self.check_type(self.configuration_id           , int                  , "configuration_id")
        self.check_type(self.finite_element_number      , int                  , "finite_element_number")
        self.check_type(self.quadrature_point_number    , int                  , "quadrature_point_number")
        self.check_type(self.quadrature_nodes           , np.ndarray           , "quadrature_nodes")
        self.check_type(self.quadrature_weights         , np.ndarray           , "quadrature_weights")
        self.check_type(self.features                   , np.ndarray           , "features")
        self.check_type(self.features_physical          , np.ndarray           , "features_physical")
        self.check_type(self.energy_density_transformed , np.ndarray           , "energy_density_transformed")
        self.check_type(self.weights_for_energy_density , np.ndarray           , "weights_for_energy_density")
        self.check_type(self.potential_transformed      , Optional[np.ndarray] , "potential_transformed")
        self.check_type(self.weights_for_potential      , Optional[np.ndarray] , "weights_for_potential")
        self.check_type(self.derivative_matrix          , Optional[np.ndarray] , "derivative_matrix")
        
        # Dimension checks
        self.check_dim(self.quadrature_nodes           , 1, "quadrature_nodes")
        self.check_dim(self.quadrature_weights         , 1, "quadrature_weights")
        self.check_dim(self.features                   , 2, "features")
        self.check_dim(self.features_physical          , 2, "features_physical")
        self.check_dim(self.energy_density_transformed , 2, "energy_density_transformed")
        self.check_dim(self.weights_for_energy_density , 2, "weights_for_energy_density")
        if self.derivative_matrix is not None:
            self.check_dim(self.derivative_matrix, 3, "derivative_matrix")


        # Shape consistency checks
        n_grid = self.quadrature_nodes.shape[0]
        n_grid_filtered = self.features.shape[0]
        assert self.finite_element_number * self.quadrature_point_number == n_grid, \
            f"finite_element_number ({self.finite_element_number}) * quadrature_point_number ({self.quadrature_point_number}) must equal n_grid ({n_grid})"
        assert self.quadrature_weights.shape[0] == n_grid, \
            f"quadrature_weights shape {self.quadrature_weights.shape} must match quadrature_nodes first dimension {n_grid}"
        assert self.features_physical.shape == self.features.shape, \
            f"features_physical shape {self.features_physical.shape} must match features shape {self.features.shape}"
        assert self.energy_density_transformed.shape[0] == n_grid_filtered, \
            f"energy_density_transformed shape {self.energy_density_transformed.shape} must match features first dimension (n_grid_filtered={n_grid_filtered})"
        assert self.weights_for_energy_density.shape[0] == n_grid_filtered, \
            f"weights_for_energy_density shape {self.weights_for_energy_density.shape} must match features first dimension (n_grid_filtered={n_grid_filtered})"
        
        if self.potential_transformed is not None:
            self.check_dim(self.potential_transformed, 2, "potential_transformed")
            assert self.potential_transformed.shape[0] == n_grid_filtered, \
                f"potential_transformed shape {self.potential_transformed.shape} must match features first dimension (n_grid_filtered={n_grid_filtered})"
        
        if self.weights_for_potential is not None:
            self.check_dim(self.weights_for_potential, 2, "weights_for_potential")
            assert self.weights_for_potential.shape[0] == n_grid_filtered, \
                f"weights_for_potential shape {self.weights_for_potential.shape} must match features first dimension (n_grid_filtered={n_grid_filtered})"
        


    def print_info(self, label: str = "Exchange-Correlation Energy Configuration"):
        """
        Print information about this ExcConfiguration.
        configuration_id: 0 = main (converged), 1..N = outer_iter_N, -1 = unknown (see data_loading config_id).
        """
        # Config type from configuration_id
        if self.configuration_id == 0:
            config_type = "main (converged)"
        elif self.configuration_id > 0:
            config_type = "outer_iter_{}".format(self.configuration_id)
        else:
            config_type = "unknown (config_id=-1)"
        
        # Derivative matrix source
        derivative_source = "stored in config" if self.derivative_matrix is not None else "use shared from loader"

        print(f"\n{'='*75}")
        print(f"{label} Summary".center(75))
        print(f"{'='*75}")
        print(f"\t Index                      : {self.index}")
        print(f"\t n_grid                     : {self.n_grid}")
        print(f"\t n_grid_filtered            : {self.n_grid_filtered}")
        print(f"\t max_element_cutoff_index   : {self.max_element_cutoff_index}")
        print(f"\t finite_element_number      : {self.finite_element_number}")
        print(f"\t quadrature_point_number    : {self.quadrature_point_number}")
        print(f"\t total_energy_true          : {self.total_energy_true if self.total_energy_true is not None else 'None'}")
        print(f"\t Config type                : {config_type}")
        print(f"\t Derivative matrix source   : {derivative_source}")
        print()
        print(f"\t quadrature_nodes           : Array of shape {self.quadrature_nodes.shape}")
        print(f"\t quadrature_weights         : Array of shape {self.quadrature_weights.shape}")
        print(f"\t features                   : Array of shape {self.features.shape}")
        print(f"\t features_physical          : Array of shape {self.features_physical.shape}")
        print(f"\t energy_density_transformed : Array of shape {self.energy_density_transformed.shape}")
        print(f"\t weights_for_energy_density : Array of shape {self.weights_for_energy_density.shape}")
        print(f"\t potential_transformed      : {self.potential_transformed.shape if self.potential_transformed is not None else 'None'}")
        print(f"\t weights_for_potential      : {self.weights_for_potential.shape if self.weights_for_potential is not None else 'None'}")
        print(f"\t derivative_matrix          : {self.derivative_matrix.shape if self.derivative_matrix is not None else 'None'}")
        print(f"{'='*75}")


    @property
    def n_grid(self) -> int:
        return self.quadrature_nodes.shape[0]

    @property
    def n_grid_filtered(self) -> int:
        return self.features.shape[0]
    
    @property
    def max_element_cutoff_index(self) -> int:
        return int(self.n_grid_filtered / self.quadrature_point_number)
    
    @property
    def quadrature_nodes_filtered(self) -> np.ndarray:
        return self.quadrature_nodes[:self.n_grid_filtered]
    
    @property
    def quadrature_weights_filtered(self) -> np.ndarray:
        return self.quadrature_weights[:self.n_grid_filtered]



@dataclass
class ComputationConfig(TypeCheckMixin):
    """
    Configuration dataclass for energy and potential computation.
    
    This dataclass contains all the necessary parameters for computing energy density
    and XC potential from model predictions (e.g. TorchXCModel.compute_energy_density_from_features,
    TorchXCModel.compute_potential_from_energy_density).
    """
    features_list                  : List[str] = field(default_factory=list)
    target_component_for_potential : Optional[PotentialTargetComponent] = None
    target_component_for_energy    : Optional[EnergyTargetComponent] = None
    
    # optional parameters for scaling
    scale_features          : bool = True
    scale_potential         : bool = True
    scale_energy            : bool = True
    scaler_for_features     : Optional[ScalerType] = None
    scaler_for_potential    : Optional[ScalerType] = None
    scaler_for_energy       : Optional[ScalerType] = None
    
    # optional parameters for symlog transformation
    use_symlog_for_features : bool = True
    use_symlog_for_potential: bool = True
    use_symlog_for_energy   : bool = True
    linthresh_for_features  : Optional[float] = 0.002
    linthresh_for_potential : Optional[float] = 0.002
    linthresh_for_energy    : Optional[float] = 0.002


    def __post_init__(self):
        # type checks
        self.check_type(self.features_list            , List[str]            , "features_list")
        self.check_type(self.scale_features           , bool                 , "scale_features")
        self.check_type(self.scale_potential          , bool                 , "scale_potential")
        self.check_type(self.scale_energy             , bool                 , "scale_energy")
        self.check_type(self.scaler_for_features      , Optional[ScalerType] , "scaler_for_features")
        self.check_type(self.scaler_for_potential     , Optional[ScalerType] , "scaler_for_potential")
        self.check_type(self.scaler_for_energy        , Optional[ScalerType] , "scaler_for_energy")

        # check if "rho" is in features_list
        assert "rho" in self.features_list, \
            FEATURES_LIST_MUST_CONTAIN_RHO_ERROR.format(self.features_list)

    @property
    def rho_index(self) -> int:
        return self.features_list.index("rho")

    def to_dict(self, exclude_scalers: bool = False) -> Dict[str, Any]:
        """Convert to dict for JSON serialization or model init. Set exclude_scalers=True for JSON."""
        from dataclasses import asdict
        d = asdict(self)
        if exclude_scalers:
            for key in ("scaler_for_features", "scaler_for_potential", "scaler_for_energy"):
                if key in d:
                    d[key] = None
        return d


@dataclass
class ExcDataLoader(TypeCheckMixin):
    """
    Data class for E_xc data loader. Data is organized by **n_configurations** (one entry per
    atomic configuration), not by n_samples (flattened grid points). Each configuration stores
    Energy (required) and optionally Potential (per-grid) data.
    """

    include_potential                : bool # Whether potential (potential_transformed, etc.) is included; prefer True when available
    atomic_numbers_per_configuration : np.ndarray  # (n_configurations,) 
    n_electrons_per_configuration    : np.ndarray  # (n_configurations,) 
    configuration_data_list          : List[ExcConfiguration] = field(default_factory=list)
    shared_derivative_matrix         : Optional[np.ndarray] = None  # Shared derivative matrix for configurations where derivative_matrix is None

    # parameters for documentation
    target_functional                : str = ""
    target_component                 : EnergyTargetComponent = ""
    target_component_for_potential   : Optional[str] = None  # When set (e.g. "v_xc"), used for potential target/config instead of target_component
    target_mode                      : str = "absolute"
    reference_functional             : Optional[str] = None
    features_list                    : List[str] = field(default_factory=list)

    # optional parameters for scaling
    scale_features                   : bool = True
    scale_potential                  : bool = True
    scale_energy                     : bool = True
    scaler_type_for_features         : str = 'robust'
    scaler_type_for_potential        : str = 'robust'
    scaler_type_for_energy           : str = 'robust'
    scaler_kwargs_for_features       : Dict[str, Any] = field(default_factory=dict)
    scaler_kwargs_for_potential      : Dict[str, Any] = field(default_factory=dict)
    scaler_kwargs_for_energy         : Dict[str, Any] = field(default_factory=dict)
    scaler_for_features              : Optional[ScalerType] = None
    scaler_for_potential             : Optional[ScalerType] = None
    scaler_for_energy                : Optional[ScalerType] = None

    # optional parameters for symlog transformation
    use_symlog_for_features          : bool = True
    use_symlog_for_potential         : bool = True
    use_symlog_for_energy            : bool = True
    linthresh_for_features           : Optional[float] = 0.002
    linthresh_for_potential          : Optional[float] = 0.002
    linthresh_for_energy             : Optional[float] = 0.002

    loss_norm_for_potential  : LossNorm = "L1norm"
    loss_norm_for_energy     : LossNorm = "L1norm"

    def __post_init__(self):
        # -------------------------------------------------------------------------
        # Type checks
        # -------------------------------------------------------------------------
        self.check_type(self.configuration_data_list          , list                , "configuration_data_list")
        self.check_type(self.include_potential                , bool                , "include_potential")
        self.check_type(self.atomic_numbers_per_configuration , np.ndarray          , "atomic_numbers_per_configuration")
        self.check_type(self.n_electrons_per_configuration    , np.ndarray          , "n_electrons_per_configuration")
        self.check_type(self.target_functional                , str                 , "target_functional")
        self.check_type(self.target_component                 , str                 , "target_component")
        self.check_type(self.target_mode                      , str                 , "target_mode")
        self.check_type(self.reference_functional             , Optional[str]       , "reference_functional")
        self.check_type(self.features_list                    , List[str]           , "features_list")
        self.check_type(self.scale_features                   , bool                , "scale_features")
        self.check_type(self.scale_potential                  , bool                , "scale_potential")
        self.check_type(self.scale_energy                     , bool                , "scale_energy")
        self.check_type(self.scaler_type_for_features         , str                 , "scaler_type_for_features")
        self.check_type(self.scaler_type_for_potential        , str                 , "scaler_type_for_potential")
        self.check_type(self.scaler_type_for_energy           , str                 , "scaler_type_for_energy")
        self.check_type(self.scaler_kwargs_for_features       , Dict[str, Any]      , "scaler_kwargs_for_features")
        self.check_type(self.scaler_kwargs_for_potential      , Dict[str, Any]      , "scaler_kwargs_for_potential")
        self.check_type(self.scaler_kwargs_for_energy         , Dict[str, Any]      , "scaler_kwargs_for_energy")
        self.check_type(self.scaler_for_features              , ScalerType          , "scaler_for_features")
        self.check_type(self.scaler_for_potential             , ScalerType          , "scaler_for_potential")
        self.check_type(self.scaler_for_energy                , ScalerType          , "scaler_for_energy")
        self.check_type(self.use_symlog_for_features          , bool                , "use_symlog_for_features")
        self.check_type(self.use_symlog_for_potential         , bool                , "use_symlog_for_potential")
        self.check_type(self.use_symlog_for_energy            , bool                , "use_symlog_for_energy")
        self.check_type(self.linthresh_for_features           , Optional[float]     , "linthresh_for_features")
        self.check_type(self.linthresh_for_potential          , Optional[float]     , "linthresh_for_potential")
        self.check_type(self.linthresh_for_energy             , Optional[float]     , "linthresh_for_energy")
        self.check_type(self.shared_derivative_matrix         , Optional[np.ndarray], "shared_derivative_matrix")
        self.check_type(self.loss_norm_for_potential          , str                 , "loss_norm_for_potential")
        self.check_type(self.loss_norm_for_energy             , str                 , "loss_norm_for_energy")

        # -------------------------------------------------------------------------
        # Configuration list structure: each element must be an ExcConfiguration instance
        # -------------------------------------------------------------------------
        assert all(isinstance(cfg, ExcConfiguration) for cfg in self.configuration_data_list), \
            "configuration_data_list elements must be ExcConfiguration instances."
        if self.include_potential:
            assert all(cfg.potential_transformed is not None for cfg in self.configuration_data_list), \
                "include_potential=True requires 'potential_transformed' in every configuration."
            assert all(cfg.weights_for_potential is not None for cfg in self.configuration_data_list), \
                "include_potential=True requires 'weights_for_potential' in every configuration."

        # -------------------------------------------------------------------------
        # Dimension and shape checks: per-configuration arrays (length = n_configurations)
        # -------------------------------------------------------------------------
        self.check_dim(self.atomic_numbers_per_configuration, 1, "atomic_numbers_per_configuration")
        self.check_dim(self.n_electrons_per_configuration   , 1, "n_electrons_per_configuration")
        n_cfg = len(self.configuration_data_list)
        self.check_shape(self.atomic_numbers_per_configuration, (n_cfg,), "atomic_numbers_per_configuration")
        self.check_shape(self.n_electrons_per_configuration   , (n_cfg,), "n_electrons_per_configuration")

        # -------------------------------------------------------------------------
        # Value checks
        # -------------------------------------------------------------------------
        assert self.target_mode in ["absolute", "delta"], \
            "target_mode must be 'absolute' or 'delta'."
        if self.scale_features:
            self.check_is_not_none(self.scaler_for_features, "scaler_for_features", "scale_features", True)
        if self.scale_potential:
            self.check_is_not_none(self.scaler_for_potential, "scaler_for_potential", "scale_potential", True)
        if self.scale_energy and self.scaler_for_energy is not None:
            self.check_is_not_none(self.scaler_for_energy, "scaler_for_energy", "scale_energy", True)
        if self.use_symlog_for_features:
            self.check_is_not_none(self.linthresh_for_features, "linthresh_for_features", "use_symlog_for_features", True)
        if self.use_symlog_for_potential:
            self.check_is_not_none(self.linthresh_for_potential, "linthresh_for_potential", "use_symlog_for_potential", True)
        if self.use_symlog_for_energy:
            self.check_is_not_none(self.linthresh_for_energy, "linthresh_for_energy", "use_symlog_for_energy", True)
        _validate_loss_norm("loss_norm_for_potential", self.loss_norm_for_potential)
        _validate_loss_norm("loss_norm_for_energy", self.loss_norm_for_energy)



    # ---------- Data loader methods (Needed for torch dataLoader) ----------
    def __len__(self) -> int:
        return self.n_configurations
    

    def __getitem__(self, index: int) -> Dict[str, Any]:
        """
        Return a configuration dict with numpy arrays (will be converted to torch tensors by collate_fn).
        """
        return self.configuration_data_list[index]


    # ---------- Primary access: by configuration ----------
    def get_configuration(self, index: int) -> Dict[str, Any]:
        """Return the configuration dict at index."""
        return self.configuration_data_list[index]


    def print_info(self, label: str = "Exc Data Loader"):
        print(f"\n{'='*75}")
        print(f"{label} Summary".center(75))
        print(f"{'='*75}")
        print(f"\tNumber of atoms             : {self.n_atoms}")
        print(f"\tNumber of configurations    : {self.n_configurations}")
        print(f"\tTotal samples (grid points) : {self.n_samples}")
        print(f"\tInclude potential           : {self.include_potential}")
        print(f"\tTarget functional           : {self.target_functional}")
        print(f"\tTarget component            : {self.target_component}")
        print(f"\tTarget mode                 : {self.target_mode}")
        print(f"\tReference functional        : {self.reference_functional}")
        print(f"\tFeatures list               : {self.features_list}")

        def _fmt_sci(x):
            a = np.asarray(x)
            return np.array2string(a, formatter={"float_kind": lambda v: f"{v: .3e}"})

        print(f"\tUse symlog for features     : {self.use_symlog_for_features}")
        print(f"\tUse symlog for potential    : {self.use_symlog_for_potential}")
        print(f"\tUse symlog for energy       : {self.use_symlog_for_energy}")
        print(f"\tLinthresh for features      : {self.linthresh_for_features}")
        print(f"\tLinthresh for potential     : {self.linthresh_for_potential}")
        print(f"\tLinthresh for energy        : {self.linthresh_for_energy}")
        print(f"\tloss_norm_for_energy        : {self.loss_norm_for_energy}")
        if self.include_potential:
            print(f"\tloss_norm_for_potential     : {self.loss_norm_for_potential}")
        if self.scale_features:
            print(f"\tScale features              : {self.scale_features}")
            print(f"\t    - mean = {_fmt_sci(self.scaler_for_features.center_)}")
            print(f"\t    - std  = {_fmt_sci(self.scaler_for_features.scale_)}")

        if self.scale_potential:
            print(f"\tScale potential             : {self.scale_potential}")
            print(f"\t    - mean = {_fmt_sci(self.scaler_for_potential.center_)}")
            print(f"\t    - std  = {_fmt_sci(self.scaler_for_potential.scale_)}")
        
        if self.scale_energy:
            print(f"\tScale energy                : {self.scale_energy}")
            print(f"\t    - mean = {_fmt_sci(self.scaler_for_energy.center_)}")
            print(f"\t    - std  = {_fmt_sci(self.scaler_for_energy.scale_)}")
        
        print(f"{'='*75}")


    @property
    def n_configurations(self) -> int:
        return len(self.configuration_data_list)


    @property
    def n_samples(self) -> int:
        return sum(cfg.features.shape[0] for cfg in self.configuration_data_list)



    @property
    def configuration_ids_per_sample(self) -> Optional[np.ndarray]:
        parts = []
        for cfg in self.configuration_data_list:
            n = cfg.features.shape[0]
            cid = int(cfg.configuration_id)
            parts.append(np.full(n, cid))
        return np.concatenate(parts, axis=0) if parts else None


    @property
    def n_atoms(self) -> int:
        return len(np.unique(self.atomic_numbers_per_configuration))

    @property
    def n_features(self) -> int:
        return len(self.features_list)

    @property
    def n_targets(self) -> int:
        if not self.configuration_data_list:
            return 0
        y = self.configuration_data_list[0].energy_density_transformed
        return y.shape[1] if y.ndim > 1 else 1

    @property
    def atomic_number_list(self) -> List[int]:
        return np.unique(self.atomic_numbers_per_configuration).tolist()


    def get_computation_config(self) -> "ComputationConfig":
        """
        Get configuration for energy and potential computation.
        
        Returns a ComputationConfig dataclass containing all necessary parameters
        for computing energy density and XC potential from model predictions.
        
        Returns
        -------
        ComputationConfig
            Configuration object with all parameters needed for computation.
        """
        return ComputationConfig(
            features_list                  = self.features_list,
            target_component_for_potential = self.target_component_for_potential,
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

    def split_data_by_atom(
        self,
        test_size    : float = 0.2,
        val_size     : float = 0.1,
        random_state : int   = 42,
        ensure_train_atoms: Optional[List[int]] = None,
    ) -> Tuple["ExcDataLoader", "ExcDataLoader", "ExcDataLoader"]:

        """Split by atom (atomic number) so that train/val/test do not share atoms."""
        offset = 0
        an_per_cfg = []
        for cfg in self.configuration_data_list:
            n = cfg.features.shape[0]
            an = self.atomic_numbers_per_configuration[offset:offset + n]
            an_per_cfg.append(int(np.unique(an)[0]))
            offset += n
        an_per_cfg = np.array(an_per_cfg)
        unique_atoms = np.unique(an_per_cfg)
        if ensure_train_atoms is not None:
            ensure_train_atoms = np.asarray(ensure_train_atoms)
            atoms_guaranteed_train = np.intersect1d(unique_atoms, ensure_train_atoms)
            atoms_to_split = np.setdiff1d(unique_atoms, atoms_guaranteed_train)
        else:
            atoms_guaranteed_train = np.array([], dtype=int)
            atoms_to_split = unique_atoms
        if len(atoms_to_split) > 0:
            atoms_train_temp, atoms_temp = train_test_split(
                atoms_to_split, test_size=(test_size + val_size), random_state=random_state
            )
            val_ratio = val_size / (test_size + val_size)
            atoms_val, atoms_test = train_test_split(
                atoms_temp, test_size=(1 - val_ratio), random_state=random_state
            )
            atoms_train = np.concatenate([atoms_guaranteed_train, atoms_train_temp])
        else:
            atoms_train = atoms_guaranteed_train
            atoms_val   = np.array([], dtype=int)
            atoms_test  = np.array([], dtype=int)
        train_mask = np.isin(an_per_cfg, atoms_train)
        val_mask   = np.isin(an_per_cfg, atoms_val)
        test_mask  = np.isin(an_per_cfg, atoms_test)

        def _subset(mask: np.ndarray) -> "ExcDataLoader":
            configs = [c for c, m in zip(self.configuration_data_list, mask) if m]
            # Per-configuration arrays for the subset (length = len(configs))
            an_per_cfg = self.atomic_numbers_per_configuration[mask]
            ne_per_cfg = self.n_electrons_per_configuration[mask]
            return ExcDataLoader(
                configuration_data_list          = configs,
                include_potential                = self.include_potential,
                atomic_numbers_per_configuration = an_per_cfg,
                n_electrons_per_configuration    = ne_per_cfg,
                shared_derivative_matrix         = self.shared_derivative_matrix,
                target_functional                = self.target_functional,
                target_component                 = self.target_component,
                target_component_for_potential   = self.target_component_for_potential,
                target_mode                      = self.target_mode,
                reference_functional             = self.reference_functional,
                features_list                    = self.features_list,
                scale_features                   = self.scale_features,
                scale_potential                  = self.scale_potential,
                scale_energy                     = self.scale_energy,
                scaler_type_for_features         = self.scaler_type_for_features,
                scaler_type_for_potential        = self.scaler_type_for_potential,
                scaler_type_for_energy           = self.scaler_type_for_energy,
                scaler_kwargs_for_features       = self.scaler_kwargs_for_features,
                scaler_kwargs_for_potential      = self.scaler_kwargs_for_potential,
                scaler_kwargs_for_energy         = self.scaler_kwargs_for_energy,
                scaler_for_features              = self.scaler_for_features,
                scaler_for_potential             = self.scaler_for_potential,
                scaler_for_energy                = self.scaler_for_energy,
                use_symlog_for_features          = self.use_symlog_for_features,
                use_symlog_for_potential         = self.use_symlog_for_potential,
                use_symlog_for_energy            = self.use_symlog_for_energy,
                linthresh_for_features           = self.linthresh_for_features,
                linthresh_for_potential          = self.linthresh_for_potential,
                linthresh_for_energy             = self.linthresh_for_energy,
                loss_norm_for_potential          = self.loss_norm_for_potential,
                loss_norm_for_energy             = self.loss_norm_for_energy,
            )

        return _subset(train_mask), _subset(val_mask), _subset(test_mask)



@dataclass
class AtomicDataset:
    """
    Data class for atomic dataset.

    Parameters
    ----------
    data_root                       : str, Dataset root directory
    scf_xc_functional               : str, SCF XC functional
    forward_pass_xc_functional_list : List[str], Forward pass XC functional list
    features_list                   : List[str], Features list
    radius_cutoff_rho_threshold     : float, Radius cutoff rho threshold
    radius_cutoff_v_x_threshold     : float, Radius cutoff v_x threshold
    radius_cutoff_v_c_threshold     : float, Radius cutoff v_c threshold
    smooth_radius_threshold         : Optional[float], Smooth radius threshold
    smooth_method                   : Optional[str], Smooth method
    smooth_kwargs                   : Dict[str, Any], Smooth kwargs
    configuration_data_list         : List[SingleConfigurationData], Configuration data list
    shared_derivative_matrix        : Optional[np.ndarray], Shared derivative matrix stored at dataset root
    """

    # Basic attributes
    data_root                       : str
    scf_xc_functional               : str
    forward_pass_xc_functional_list : List[str]
    features_list                   : List[str]

    # Other attributes
    radius_cutoff_rho_threshold     : float           = 1e-6
    radius_cutoff_v_x_threshold     : float           = 1e-8
    radius_cutoff_v_c_threshold     : float           = 1e-8
    smooth_radius_threshold         : Optional[float] = None
    smooth_method                   : Optional[str]   = None
    smooth_kwargs                   : Dict[str, Any]  = field(default_factory=dict)

    # Data attributes
    configuration_data_list         : List[SingleConfigurationData] = field(default_factory=list)
    shared_derivative_matrix        : Optional[np.ndarray] = None  # Shared derivative matrix stored at dataset root


    def __post_init__(self):

        # Type and value checks
        AtomicDataManager.check_data_root(self.data_root)
        AtomicDataManager.check_features_list(self.features_list)
        AtomicDataManager.check_scf_xc_functional(self.scf_xc_functional)
        AtomicDataManager.check_forward_pass_xc_functional_list(self.forward_pass_xc_functional_list)
        
        AtomicDataManager.check_radius_cutoff_threshold(
            self.radius_cutoff_rho_threshold,
            self.radius_cutoff_v_x_threshold,
            self.radius_cutoff_v_c_threshold,
        )
        AtomicDataManager.check_configuration_data_list(self.configuration_data_list)


        # Initialize data attributes
        self.v_xc_is_smoothed                  = False
        self.e_xc_is_smoothed                  = False
        self.potential_weights_data_is_updated = False
        self.energy_weights_data_is_updated    = False
        self.include_energy_density = self._get_include_energy_density()

        # Cached data attributes
        self.cached_atomic_numbers_unique            : Optional[np.ndarray] = None
        self.cached_atomic_numbers_per_configuration : Optional[np.ndarray] = None
        self.cached_atomic_numbers_per_sample        : Optional[np.ndarray] = None
        self.cached_n_electrons_per_sample           : Optional[np.ndarray] = None
    
        self.cached_cutoff_radii                     : Optional[np.ndarray] = None
        self.cached_cutoff_indices                   : Optional[np.ndarray] = None
        self.cached_quadrature_nodes                 : Optional[np.ndarray] = None
        self.cached_configuration_ids                : Optional[List[int]]  = None

        self.cached_features_data                    : Optional[np.ndarray] = None
        self.cached_scf_xc_data                      : Optional[XCDataType] = None
        self.cached_forward_pass_xc_data_list        : Optional[List[XCDataType]] = None


    def _get_include_energy_density(self) -> bool:
        """
        Check if energy density is included.
        """
        include_energy_density = True

        # Check if energy density is included in SCF XC data
        for configuration_data in self.configuration_data_list:
            if configuration_data.scf_xc_data[2] is None or configuration_data.scf_xc_data[3] is None:
                include_energy_density = False
                break
        
        # Check if energy density is included in forward pass XC data
        for configuration_data in self.configuration_data_list:
            for forward_pass_xc_data in configuration_data.forward_pass_xc_data_list:
                if forward_pass_xc_data[2] is None or forward_pass_xc_data[3] is None:
                    include_energy_density = False
                    break

        return include_energy_density


    def smooth_xc_data(
        self,
        smooth_radius_threshold : float = 10.0,
        smooth_method           : str   = 'savgol',
        smooth_kwargs           : Dict[str, Any] = field(default_factory=dict),
    ) -> None:
        """
        Smooth V_xc (potential) and E_xc (energy density) data for r > r_threshold.
        """
        # Type and value check
        AtomicDataManager.check_smooth_parameters(smooth_radius_threshold, smooth_method, smooth_kwargs)

        # Update smoothing parameters (shared by both V_xc and E_xc)
        self.smooth_radius_threshold = smooth_radius_threshold
        self.smooth_method           = smooth_method
        self.smooth_kwargs           = smooth_kwargs

        # Helper function to smooth potential (v_x, v_c)
        def smooth_potential(v_x, v_c, quadrature_nodes):
            v_x = DataProcessor.smooth_radial_data(v_x, quadrature_nodes, self.smooth_radius_threshold, self.smooth_method, **self.smooth_kwargs)
            v_c = DataProcessor.smooth_radial_data(v_c, quadrature_nodes, self.smooth_radius_threshold, self.smooth_method, **self.smooth_kwargs)
            return v_x, v_c

        # Helper function to smooth energy density (e_x, e_c)
        def smooth_energy(e_x, e_c, quadrature_nodes):
            e_x = DataProcessor.smooth_radial_data(e_x, quadrature_nodes, self.smooth_radius_threshold, self.smooth_method, **self.smooth_kwargs)
            e_c = DataProcessor.smooth_radial_data(e_c, quadrature_nodes, self.smooth_radius_threshold, self.smooth_method, **self.smooth_kwargs)
            return e_x, e_c

        smoothed_this_call = False

        # Smooth V_xc data if not already smoothed
        if not self.v_xc_is_smoothed:
            for configuration_data in self.configuration_data_list:
                # Smooth SCF XC potential data
                scf_v_x, scf_v_c = smooth_potential(configuration_data.scf_xc_data[0], configuration_data.scf_xc_data[1], configuration_data.quadrature_nodes)
                configuration_data.scf_xc_data = (scf_v_x, scf_v_c, configuration_data.scf_xc_data[2], configuration_data.scf_xc_data[3])

                # Smooth forward pass XC potential data
                for idx, forward_pass_xc_data in enumerate(configuration_data.forward_pass_xc_data_list):
                    fp_v_x, fp_v_c = smooth_potential(forward_pass_xc_data[0], forward_pass_xc_data[1], configuration_data.quadrature_nodes)
                    configuration_data.forward_pass_xc_data_list[idx] = (fp_v_x, fp_v_c, forward_pass_xc_data[2], forward_pass_xc_data[3])

            self.v_xc_is_smoothed = True
            smoothed_this_call = True
        else:
            print(V_XC_IS_ALREADY_SMOOTHED_WARNING)

        # Smooth E_xc data if not already smoothed and energy density is included
        if self.include_energy_density and not self.e_xc_is_smoothed:
            for configuration_data in self.configuration_data_list:
                e_x, e_c = configuration_data.scf_xc_data[2], configuration_data.scf_xc_data[3]
                if e_x is not None and e_c is not None:
                    scf_e_x, scf_e_c = smooth_energy(e_x, e_c, configuration_data.quadrature_nodes)
                    configuration_data.scf_xc_data = (configuration_data.scf_xc_data[0], configuration_data.scf_xc_data[1], scf_e_x, scf_e_c)

                for idx, forward_pass_xc_data in enumerate(configuration_data.forward_pass_xc_data_list):
                    fp_e_x, fp_e_c = forward_pass_xc_data[2], forward_pass_xc_data[3]
                    if fp_e_x is not None and fp_e_c is not None:
                        fp_e_x_smooth, fp_e_c_smooth = smooth_energy(fp_e_x, fp_e_c, configuration_data.quadrature_nodes)
                        v_x, v_c = configuration_data.forward_pass_xc_data_list[idx][0], configuration_data.forward_pass_xc_data_list[idx][1]
                        configuration_data.forward_pass_xc_data_list[idx] = (v_x, v_c, fp_e_x_smooth, fp_e_c_smooth)

            self.e_xc_is_smoothed = True
            smoothed_this_call = True
        elif self.include_energy_density:
            print(E_XC_IS_ALREADY_SMOOTHED_WARNING)

        # Clean the cached data and weights if any smoothing was done
        if smoothed_this_call:
            self.cached_scf_xc_data = None
            self.cached_forward_pass_xc_data_list = None

            if self.potential_weights_data_is_updated or self.energy_weights_data_is_updated:
                for configuration_data in self.configuration_data_list:
                    configuration_data.clear_weights_data()
                self.potential_weights_data_is_updated = False
                self.energy_weights_data_is_updated    = False


    def clip_positive_energy_density(self, epsilon: float = 1E-35) -> None:
        """
        Set all positive XC energy density values to 0 (in-place).
        
        Due to numerical errors, some exchange-correlation energy densities (e_x, e_c)
        may be slightly positive when they should be negative. This method clips
        all values > 0 to 0 for both SCF and forward-pass XC data.
        
        Requires include_energy_density=True. Clears cached XC data and weights.
        """
        if not self.include_energy_density:
            raise ValueError(
                "clip_positive_energy_density requires include_energy_density=True. "
                "The dataset does not contain energy density data."
            )
        
        for configuration_data in self.configuration_data_list:
            # Clip SCF energy density (e_x, e_c)
            e_x, e_c = configuration_data.scf_xc_data[2], configuration_data.scf_xc_data[3]
            if e_x is not None:
                e_x[e_x > -epsilon] = -epsilon
            if e_c is not None:
                e_c[e_c > -epsilon] = -epsilon
            
            # Clip forward-pass energy density for each functional
            for idx, forward_pass_xc_data in enumerate(configuration_data.forward_pass_xc_data_list):
                fp_e_x, fp_e_c = forward_pass_xc_data[2], forward_pass_xc_data[3]
                if fp_e_x is not None:
                    fp_e_x[fp_e_x > -epsilon] = -epsilon
                if fp_e_c is not None:
                    fp_e_c[fp_e_c > -epsilon] = -epsilon
        
        # Clear cached data
        self.cached_scf_xc_data = None
        self.cached_forward_pass_xc_data_list = None
        
        if self.potential_weights_data_is_updated or self.energy_weights_data_is_updated:
            for configuration_data in self.configuration_data_list:
                configuration_data.clear_weights_data()
            self.potential_weights_data_is_updated = False
            self.energy_weights_data_is_updated    = False


    def print_info(self):
        """
        Print information about the atomic dataset.
        """

        print(f"{'='*75}")
        print(f"Loaded Data Summary".center(75))
        print(f"{'='*75}")
        
        print(f"\tNumber of atoms                    : {self.n_atoms}")
        print(f"\tNumber of configurations           : {self.n_configurations}")
        print(f"\tNumber of samples                  : {self.n_samples}")
        print(f"\tsmoothed XC data                   : {self.v_xc_is_smoothed and self.e_xc_is_smoothed}")
        print(f"\tInclude energy density             : {self.include_energy_density}")
        print(f"\tCutoff radii range                 : [{np.min(self.cutoff_radii):.6f}, {np.max(self.cutoff_radii):.6f}]")
        print(f"\tSCF XC functional                  : {self.scf_xc_functional}")
        print(f"\tForward pass XC functional list    : {self.forward_pass_xc_functional_list}")
        print(f"\tFeatures in features list          : {len(self.features_list)} channels")
        for idx, feature in enumerate(self.features_list):
            suffix = " (repeated)" if feature in self.features_list[:idx] else ""
            print(f"\t    - Channel {idx + 1}: {feature}{suffix}")

        print()
        print(f"\tradius_cutoff_rho_threshold        : {self.radius_cutoff_rho_threshold}")
        print(f"\tradius_cutoff_v_x_threshold        : {self.radius_cutoff_v_x_threshold}")
        print(f"\tradius_cutoff_v_c_threshold        : {self.radius_cutoff_v_c_threshold}")
        print(f"\tv_xc is smoothed                   : {self.v_xc_is_smoothed}")
        print(f"\te_xc is smoothed                   : {self.e_xc_is_smoothed}")
        print(f"\tsmooth_radius_threshold            : {self.smooth_radius_threshold}")
        print(f"\tsmooth_method                      : {self.smooth_method}")
        print(f"\tsmooth_kwargs                      : {self.smooth_kwargs}")
        print(f"\tpotential_weights_data_is_updated  : {self.potential_weights_data_is_updated}")
        print(f"\tenergy_weights_data_is_updated     : {self.energy_weights_data_is_updated}")
        print(f"\tinclude_energy_density             : {self.include_energy_density}")
        
        # Print derivative matrix information
        print(f"\tShared derivative matrix exists    : {self.shared_derivative_matrix is not None}")
        if self.shared_derivative_matrix is not None:
            print(f"\tShared derivative matrix shape     : {self.shared_derivative_matrix.shape}")

        print()
        print(f"\tshape of cutoff_radii              : Array of shape {self.cutoff_radii.shape}")
        print(f"\tshape of cutoff_indices            : Array of shape {self.cutoff_indices.shape}")
        print(f"\tshape of atomic_numbers_per_sample : Array of shape {self.atomic_numbers_per_sample.shape}")
        print(f"\tshape of quadrature_nodes          : Array of shape {self.quadrature_nodes.shape}")
        print(f"\tshape of configuration_ids         : List of {len(self.configuration_ids)} integers")
        print(f"\tshape of features_data             : {self.features_data.shape}")
        print(f"\tshape of scf_xc_data               : Tuple of 4 elements")
        print(f"\t    - v_x: Array of shape {self.scf_xc_data[0].shape}")
        print(f"\t    - v_c: Array of shape {self.scf_xc_data[1].shape}")
        print( "\t    - e_x: {}".format(f"Array of shape {self.scf_xc_data[2].shape}" if self.scf_xc_data[2] is not None else "None"))
        print( "\t    - e_c: {}".format(f"Array of shape {self.scf_xc_data[3].shape}" if self.scf_xc_data[3] is not None else "None"))
        print(f"\tshape of forward_pass_xc_data_list : List of {len(self.forward_pass_xc_data_list)} tuples, each with 4 elements")
        if len(self.forward_pass_xc_data_list) > 0:
            print(f"\t    - v_x: Array of shape {self.forward_pass_xc_data_list[0][0].shape}")
            print(f"\t    - v_c: Array of shape {self.forward_pass_xc_data_list[0][1].shape}")
            print( "\t    - e_x: {}".format(f"Array of shape {self.forward_pass_xc_data_list[0][2].shape}" if self.forward_pass_xc_data_list[0][2] is not None else "None"))
            print( "\t    - e_c: {}".format(f"Array of shape {self.forward_pass_xc_data_list[0][3].shape}" if self.forward_pass_xc_data_list[0][3] is not None else "None"))
        else:
            print("\t    - No forward pass XC data")
        

        
        print(f"{'='*75}")


    # Data preparation methods, for potential data
    def prepare_potential_dataloader(
        self,
        # required parameters
        target_functional              : str,
        target_component               : PotentialTargetComponent,
        reference_functional           : Optional[str],

        # optional parameters for scaling
        scale_features                 : bool = True,
        scale_potential                : bool = True,
        scaler_type_for_features       : str = 'robust',
        scaler_type_for_potential      : str = 'robust',
        scaler_kwargs_for_features     : Optional[Dict[str, Any]] = None,
        scaler_kwargs_for_potential    : Optional[Dict[str, Any]] = None,

        # optional parameters for symlog transformation
        use_symlog_for_features        : bool = True,
        use_symlog_for_potential       : bool = True,
        linthresh_for_features         : Optional[float] = 0.002,
        linthresh_for_potential        : Optional[float] = 0.002,

        # optional parameters for weights calculation
        min_weight_ratio_for_potential : float = 1e-2,
        normalize_weight_for_potential : bool = True,
        loss_norm_for_potential        : LossNorm = "L1norm",

        # optional: split into train/val/test by atom (default False = return single loader)
        split_train_val_test           : bool = False,
        test_size                      : float = 0.2,
        val_size                       : float = 0.1,
        random_state                   : int = 42,
        ensure_train_atoms             : Optional[List[int]] = None,

        # deprecated (use *_for_features / *_for_potential names instead)
        scale_targets                  : Optional[bool] = None,
        scaler_type_features           : Optional[str] = None,
        scaler_type_targets            : Optional[str] = None,
        scaler_kwargs_targets          : Optional[Dict[str, Any]] = None,
        use_symlog_features            : Optional[bool] = None,
        use_symlog_targets             : Optional[bool] = None,
        linthresh_features             : Optional[float] = None,
        linthresh_targets              : Optional[float] = None,
        min_weight_ratio               : Optional[float] = None,
    ) -> Union[VxcDataLoader, Tuple[VxcDataLoader, VxcDataLoader, VxcDataLoader]]:
        """
        Prepare data with symlog transformation and scaling, for potential data.

        Parameters
        ----------
        target_functional : str
            Target XC functional used to construct training targets
        target_component : str
            Target component to learn: "v_xc", "v_x", "v_c", or "v_x_v_c"
        reference_functional : str or None
            Reference XC functional for delta learning. If None, use absolute targets
        scale_features : bool
            Whether to scale features
        scale_potential : bool
            Whether to scale potential (target) values
        scaler_type_for_features : str
            Type of feature scaler: 'robust' or 'standard'
        scaler_type_for_potential : str
            Type of potential scaler: 'robust' or 'standard'
        use_symlog_for_features : bool
            Whether to apply symlog to features
        use_symlog_for_potential : bool
            Whether to apply symlog to potential (target) values
        linthresh_for_features, linthresh_for_potential : float
            Linear threshold for symlog
        min_weight_ratio_for_potential : float
            Minimum weight ratio for weights calculation
        normalize_weight_for_potential : bool
            Whether to normalize weights for potential
            if True, all configurations have the same importance
            if False, all grid points have the same importance
        loss_norm_for_potential : str
            ``L1norm`` (default) or ``L2norm``. Pair with ``L1Loss`` / ``MSELoss`` respectively in training.
        split_train_val_test : bool
            If True, split by atomic number: ensure_train_atoms all in training, remaining atoms split by ratio; return three loaders.
            If False (default), return a single loader on full data.
        test_size : float
            Proportion of atoms for test set (used when split_train_val_test=True).
        val_size : float
            Proportion of atoms for validation set (used when split_train_val_test=True).
        random_state : int
            Random seed for split (used when split_train_val_test=True).
        ensure_train_atoms : list or None
            Atomic numbers that must be in training set; remaining atoms split by test_size/val_size (used when split_train_val_test=True).

        Deprecated (will be removed in a future version)
        -------------------------------------------------
        scale_targets         : use scale_potential instead
        scaler_type_features  : use scaler_type_for_features instead
        scaler_type_targets   : use scaler_type_for_potential instead
        scaler_kwargs_targets : use scaler_kwargs_for_potential instead
        use_symlog_features   : use use_symlog_for_features instead
        use_symlog_targets    : use use_symlog_for_potential instead
        linthresh_features    : use linthresh_for_features instead
        linthresh_targets     : use linthresh_for_potential instead
        min_weight_ratio      : use min_weight_ratio_for_potential instead

        Returns
        -------
        VxcDataLoader or Tuple[VxcDataLoader, VxcDataLoader, VxcDataLoader]
            If split_train_val_test is False: single data loader.
            If split_train_val_test is True: (train_loader, val_loader, test_loader).
        """

        # Handle deprecated parameters
        if scale_targets is not None:
            warnings.warn(SCALE_TARGETS_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            scale_potential = scale_targets
        if scaler_type_features is not None:
            warnings.warn(SCALER_TYPE_FEATURES_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            scaler_type_for_features = scaler_type_features
        if scaler_type_targets is not None:
            warnings.warn(SCALER_TYPE_TARGETS_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            scaler_type_for_potential = scaler_type_targets
        if scaler_kwargs_targets is not None:
            warnings.warn(SCALER_KWARGS_TARGETS_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            scaler_kwargs_for_potential = scaler_kwargs_targets
        if use_symlog_features is not None:
            warnings.warn(USE_SYMLOG_FEATURES_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            use_symlog_for_features = use_symlog_features
        if use_symlog_targets is not None:
            warnings.warn(USE_SYMLOG_TARGETS_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            use_symlog_for_potential = use_symlog_targets
        if linthresh_features is not None:
            warnings.warn(LINTHRESH_FEATURES_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            linthresh_for_features = linthresh_features
        if linthresh_targets is not None:
            warnings.warn(LINTHRESH_TARGETS_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            linthresh_for_potential = linthresh_targets
        if min_weight_ratio is not None:
            warnings.warn(MIN_WEIGHT_RATIO_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            min_weight_ratio_for_potential = min_weight_ratio
        _validate_loss_norm("loss_norm_for_potential", loss_norm_for_potential)
        # Type and value checks
        if not SKLEARN_AVAILABLE:
            raise ImportError(SKLEARN_NOT_AVAILABLE_FOR_DATA_PREPROCESSING_ERROR)

        # Check if target functional exists in the dataset
        if not isinstance(target_functional, str):
            raise ValueError(TARGET_FUNCTIONAL_NOT_STRING_ERROR.format(target_functional))
        if not self.exists_functional(target_functional):
            raise ValueError(TARGET_FUNCTIONAL_NOT_IN_DATASET_ERROR.format(target_functional, self.scf_xc_functional, self.forward_pass_xc_functional_list))
        
        # check if target_component is a string and in the valid list
        if target_component not in VALID_POTENTIAL_TARGET_COMPONENTS:
            raise ValueError(TARGET_COMPONENT_NOT_IN_VALID_LIST_ERROR.format(VALID_POTENTIAL_TARGET_COMPONENTS, target_component))

        # Check if reference functional exists in the dataset
        if reference_functional is not None:
            if not isinstance(reference_functional, str):
                raise ValueError(REFERENCE_FUNCTIONAL_NOT_STRING_ERROR.format(reference_functional))
            if not self.exists_functional(reference_functional):
                raise ValueError(REFERENCE_FUNCTIONAL_NOT_IN_DATASET_ERROR.format(reference_functional, self.scf_xc_functional, self.forward_pass_xc_functional_list))

        if scaler_kwargs_for_features is None:
            scaler_kwargs_for_features = {}
        if scaler_kwargs_for_potential is None:
            scaler_kwargs_for_potential = {}

        # Determine target mode
        if reference_functional is None:
            target_mode = "absolute"
        else:
            target_mode = "delta"

        # Extract atomic numbers and number of electrons per sample
        atomic_numbers_per_sample = self.atomic_numbers_per_sample
        n_electrons_per_sample    = self.n_electrons_per_sample

        # Extract features
        X = self.features_data.copy()

        # Extract target component data
        target_component_data = self.extract_potential_component(self.get_xc_data(target_functional), target_component)
        if target_mode == "delta":
            reference_component_data = self.extract_potential_component(self.get_xc_data(reference_functional), target_component)
            y = target_component_data - reference_component_data
        else:
            y = target_component_data

        if target_component == "enhancement_factor_x":
            from ..xc.lda import lda_exchange_potential_generic
            lda_exchange_potential = lda_exchange_potential_generic(self.rho)
            y = y / lda_exchange_potential


        # Symlog and reshape y once (shared by split and no-split paths)
        if use_symlog_for_features:
            X = DataProcessor.symlog(X, linthresh=linthresh_for_features)
        if use_symlog_for_potential:
            y = DataProcessor.symlog(y, linthresh=linthresh_for_potential)
        if y.ndim == 1:
            y = y.reshape(-1, 1)

        if split_train_val_test:
            # Split first, then fit scalers only on training set (by atomic number)
            unique_atoms = np.unique(atomic_numbers_per_sample)

            # Ensure specified atoms are in training set
            if ensure_train_atoms is not None:
                ensure_train_atoms = np.array(ensure_train_atoms)
                # Find atoms that exist in the data and should be in training set
                atoms_guaranteed_train = np.intersect1d(unique_atoms, ensure_train_atoms)
                # Remaining atoms to split
                atoms_to_split = np.setdiff1d(unique_atoms, atoms_guaranteed_train)
                if len(atoms_guaranteed_train) > 0:
                    print(f"Ensuring atoms {atoms_guaranteed_train.tolist()} are in training set")
            else:
                atoms_guaranteed_train = np.array([], dtype=int)
                atoms_to_split = unique_atoms

            # Split remaining atoms into train/val/test
            if len(atoms_to_split) > 0:
                atoms_train_temp, atoms_temp = train_test_split(
                    atoms_to_split, test_size=(test_size + val_size), random_state=random_state
                )
                val_ratio = val_size / (test_size + val_size)
                atoms_val, atoms_test = train_test_split(
                    atoms_temp, test_size=(1 - val_ratio), random_state=random_state
                )
                # Combine guaranteed training atoms with randomly split training atoms
                atoms_train = np.concatenate([atoms_guaranteed_train, atoms_train_temp])
            else:
                # All atoms are guaranteed to be in training set
                atoms_train = atoms_guaranteed_train
                atoms_val   = np.array([], dtype=int)
                atoms_test  = np.array([], dtype=int)

            # Create masks (by atomic number) for train/val/test
            train_mask = np.isin(atomic_numbers_per_sample, atoms_train)
            val_mask   = np.isin(atomic_numbers_per_sample, atoms_val)
            test_mask  = np.isin(atomic_numbers_per_sample, atoms_test)

            # Subsets for fit/transform
            X_train, X_val, X_test = X[train_mask], X[val_mask], X[test_mask]
            y_train, y_val, y_test = y[train_mask], y[val_mask], y[test_mask]
            if X_val.shape[0] == 0 or X_test.shape[0] == 0:
                raise ValueError(SPLIT_EMPTY_VAL_OR_TEST_ERROR.format(X_val.shape[0], X_test.shape[0]))
            
            # Fit scalers on training set only; transform train/val/test
            scaler_X = None
            if scale_features:
                if scaler_type_for_features == 'robust':
                    scaler_X = RobustScaler(**scaler_kwargs_for_features)
                else:
                    scaler_X = StandardScaler(**scaler_kwargs_for_features)                
                scaler_X.fit(X_train)
                X_train = scaler_X.transform(X_train)
                X_val   = scaler_X.transform(X_val)
                X_test  = scaler_X.transform(X_test)
            
            scaler_y = None
            if scale_potential:
                if scaler_type_for_potential == 'robust':
                    scaler_y = RobustScaler(**scaler_kwargs_for_potential)
                else:
                    scaler_y = StandardScaler(**scaler_kwargs_for_potential)            
                scaler_y.fit(y_train)
                y_train = scaler_y.transform(y_train)
                y_val   = scaler_y.transform(y_val)
                y_test  = scaler_y.transform(y_test)

            # Weights use train-fit scaler_y (sigma_v from training set)
            self.update_potential_weights_data(
                normalize_weight       = normalize_weight_for_potential,
                min_weight_ratio       = min_weight_ratio_for_potential,
                target_component       = target_component,
                target_functional      = target_functional,
                reference_functional   = reference_functional,
                linthresh_for_targets  = linthresh_for_potential,
                scaler_y               = scaler_y,
                loss_norm              = loss_norm_for_potential,
                use_symlog_for_targets = use_symlog_for_potential,
            )
            weights_for_potential = self.potential_weights_data

            # Create loaders for train/val/test
            def _make_loader(X_part, y_part, w_part, an_part, ne_part):
                return VxcDataLoader(
                    features                    = X_part, 
                    potential_transformed       = y_part,
                    weights_for_potential       = w_part,
                    atomic_numbers_per_sample   = an_part,
                    n_electrons_per_sample      = ne_part,
                    target_functional           = target_functional,
                    target_component            = target_component,
                    target_mode                 = target_mode,
                    reference_functional        = reference_functional,
                    features_list               = self.features_list,
                    scale_features              = scale_features,
                    scale_potential             = scale_potential,
                    scaler_type_for_features    = scaler_type_for_features,
                    scaler_type_for_potential   = scaler_type_for_potential,
                    scaler_kwargs_for_features  = scaler_kwargs_for_features,
                    scaler_kwargs_for_potential = scaler_kwargs_for_potential,
                    scaler_for_features         = scaler_X,
                    scaler_for_potential        = scaler_y,
                    use_symlog_for_features     = use_symlog_for_features,
                    use_symlog_for_potential    = use_symlog_for_potential,
                    linthresh_for_features      = linthresh_for_features,
                    linthresh_for_potential     = linthresh_for_potential,
                    loss_norm_for_potential     = loss_norm_for_potential,
                )

            return (
                _make_loader(X_train, y_train, weights_for_potential[train_mask], atomic_numbers_per_sample[train_mask], n_electrons_per_sample[train_mask]),
                _make_loader(X_val  , y_val  , weights_for_potential[val_mask]  , atomic_numbers_per_sample[val_mask]  , n_electrons_per_sample[val_mask]),
                _make_loader(X_test , y_test , weights_for_potential[test_mask] , atomic_numbers_per_sample[test_mask] , n_electrons_per_sample[test_mask]),
            )

        # --- No split: fit on full data, return single loader ---
        scaler_X = None
        if scale_features:
            if scaler_type_for_features == 'robust':
                scaler_X = RobustScaler(**scaler_kwargs_for_features)
            else:
                scaler_X = StandardScaler(**scaler_kwargs_for_features)
            X = scaler_X.fit_transform(X)

        scaler_y = None
        if scale_potential:
            if scaler_type_for_potential == 'robust':
                scaler_y = RobustScaler(**scaler_kwargs_for_potential)
            else:
                scaler_y = StandardScaler(**scaler_kwargs_for_potential)
            y = scaler_y.fit_transform(y)

        self.update_potential_weights_data(
            normalize_weight       = normalize_weight_for_potential,
            min_weight_ratio       = min_weight_ratio_for_potential,
            target_component       = target_component,
            target_functional      = target_functional,
            reference_functional   = reference_functional,
            linthresh_for_targets  = linthresh_for_potential,
            scaler_y               = scaler_y,
            loss_norm              = loss_norm_for_potential,
            use_symlog_for_targets = use_symlog_for_potential,
        )
        weights_for_potential = self.potential_weights_data

        # Create loader for full data
        return VxcDataLoader(
            features                    = X, 
            potential_transformed       = y,
            weights_for_potential       = weights_for_potential,
            atomic_numbers_per_sample   = atomic_numbers_per_sample,
            n_electrons_per_sample      = n_electrons_per_sample,
            target_functional           = target_functional,
            target_component            = target_component,
            target_mode                 = target_mode,
            reference_functional        = reference_functional,
            features_list               = self.features_list,
            scale_features              = scale_features,
            scale_potential             = scale_potential,
            scaler_type_for_features    = scaler_type_for_features,
            scaler_type_for_potential   = scaler_type_for_potential,
            scaler_kwargs_for_features  = scaler_kwargs_for_features,
            scaler_kwargs_for_potential = scaler_kwargs_for_potential,
            scaler_for_features         = scaler_X,
            scaler_for_potential        = scaler_y,
            use_symlog_for_features     = use_symlog_for_features,
            use_symlog_for_potential    = use_symlog_for_potential,
            linthresh_for_features      = linthresh_for_features,
            linthresh_for_potential     = linthresh_for_potential,
            loss_norm_for_potential     = loss_norm_for_potential,
        )


    # Data preparation methods, for potential data
    def prepare_energy_dataloader(
        self,
        # required parameters
        target_functional              : str,
        target_component               : EnergyTargetComponent,
        reference_functional           : Optional[str],
        include_potential              : bool = True, # whether to include potential in the data

        # optional parameters for scaling
        scale_features                 : bool = True,
        scale_potential                : bool = True,
        scale_energy                   : bool = True,
        scaler_type_for_features       : str = 'robust',
        scaler_type_for_potential      : str = 'robust',
        scaler_type_for_energy         : str = 'robust',
        scaler_kwargs_for_features     : Optional[Dict[str, Any]] = None,
        scaler_kwargs_for_potential    : Optional[Dict[str, Any]] = None,
        scaler_kwargs_for_energy       : Optional[Dict[str, Any]] = None,

        # optional parameters for symlog transformation
        use_symlog_for_features        : bool = True,
        use_symlog_for_potential       : bool = True,
        use_symlog_for_energy          : bool = True,
        linthresh_for_features         : Optional[float] = 0.002,
        linthresh_for_potential        : Optional[float] = 0.002,
        linthresh_for_energy           : Optional[float] = 0.002,

        # optional parameters for weights calculation
        min_weight_ratio_for_potential : float = 1e-2,
        min_weight_ratio_for_energy    : float = 1e-2,
        normalize_weight_for_potential : bool  = True,
        normalize_weight_for_energy    : bool  = True,
        loss_norm_for_potential        : Union[LossNorm, str] = "L1norm",
        loss_norm_for_energy           : Union[LossNorm, str] = "L1norm",

        # optional: split into train/val/test by atom (default False = return single loader)
        split_train_val_test           : bool  = False,
        test_size                      : float = 0.2,
        val_size                       : float = 0.1,
        random_state                   : int   = 42,
        ensure_train_atoms             : Optional[List[int]] = None,
    ) -> Union[ExcDataLoader, Tuple[ExcDataLoader, ExcDataLoader, ExcDataLoader]]:
        """
        Prepare data with symlog transformation and scaling, for energy data.

        Parameters
        ----------
        target_functional : str
            Target XC functional used to construct training targets
        target_component : str
            Target component to learn: "e_xc", "e_x", "e_c", or "e_x_e_c"
        reference_functional : str or None
            Reference XC functional for delta learning. If None, use absolute targets
        scale_features : bool
            Whether to scale features
        scale_potential : bool
            Whether to scale potential (target) values
        scaler_type_for_features : str
            Type of feature scaler: 'robust' or 'standard'
        scaler_type_for_potential : str
            Type of potential scaler: 'robust' or 'standard'
        use_symlog_for_features : bool
            Whether to apply symlog to features
        use_symlog_for_potential : bool
            Whether to apply symlog to potential (target) values
        linthresh_for_features, linthresh_for_potential : float
            Linear threshold for symlog
        min_weight_ratio_for_potential : float
            Minimum weight ratio for weights calculation
        normalize_weight_for_potential : bool
            Whether to normalize weights for potential
            if True, all configurations have the same importance
            if False, all grid points have the same importance
        loss_norm_for_energy : str
            ``L1norm`` or ``L2norm`` for energy-density sample weights (see prepare_potential_dataloader).
        loss_norm_for_potential : str
            Same for potential branch when ``include_potential=True``.
        split_train_val_test : bool
            If True, split by atom: ensure_train_atoms all in training, remaining atoms split by ratio; return three loaders.
            If False (default), return a single loader on full data.
        test_size : float
            Proportion of atoms for test set (used when split_train_val_test=True).
        val_size : float
            Proportion of atoms for validation set (used when split_train_val_test=True).
        random_state : int
            Random seed for split (used when split_train_val_test=True).
        ensure_train_atoms : list or None
            Atomic numbers that must be in training set; remaining atoms split by test_size/val_size (used when split_train_val_test=True).

        Returns
        -------
        ExcDataLoader or Tuple[ExcDataLoader, ExcDataLoader, ExcDataLoader]
            If split_train_val_test is False: single data loader (by n_configurations).
            If split_train_val_test is True: (train_loader, val_loader, test_loader).
        """

        # Type and value checks
        if not SKLEARN_AVAILABLE:
            raise ImportError(SKLEARN_NOT_AVAILABLE_FOR_DATA_PREPROCESSING_ERROR)
        if not self.include_energy_density:
            raise ValueError(PREPARE_ENERGY_DATALOADER_REQUIRES_INCLUDE_ENERGY_DENSITY_ERROR)

        _validate_loss_norm("loss_norm_for_energy", loss_norm_for_energy)
        _validate_loss_norm("loss_norm_for_potential", loss_norm_for_potential)

        # check if target functional is a string and in the dataset
        if not isinstance(target_functional, str):
            raise ValueError(TARGET_FUNCTIONAL_NOT_STRING_ERROR.format(target_functional))
        if not self.exists_functional(target_functional):
            raise ValueError(TARGET_FUNCTIONAL_NOT_IN_DATASET_ERROR.format(target_functional, self.scf_xc_functional, self.forward_pass_xc_functional_list))

        # check if target_component is a string and in the valid list
        if target_component not in VALID_ENERGY_TARGET_COMPONENTS:
            raise ValueError(TARGET_COMPONENT_NOT_IN_VALID_LIST_ERROR.format(VALID_ENERGY_TARGET_COMPONENTS, target_component))
        
        # check if reference_functional is None or in the dataset
        if reference_functional is not None:
            if not isinstance(reference_functional, str):
                raise ValueError(REFERENCE_FUNCTIONAL_NOT_STRING_ERROR.format(reference_functional))
            if not self.exists_functional(reference_functional):
                raise ValueError(REFERENCE_FUNCTIONAL_NOT_IN_DATASET_ERROR.format(reference_functional, self.scf_xc_functional, self.forward_pass_xc_functional_list))

        if scaler_kwargs_for_features is None:
            scaler_kwargs_for_features = {}
        if scaler_kwargs_for_potential is None:
            scaler_kwargs_for_potential = {}
        if scaler_kwargs_for_energy is None:
            scaler_kwargs_for_energy = {}

        # determine target_mode
        if reference_functional is not None:
            target_mode = "delta"
        else:
            target_mode = "absolute"

        # Extract atomic numbers and number of electrons per configuration
        atomic_numbers_per_configuration = np.array([configuration_data.atomic_number for configuration_data in self.configuration_data_list])
        n_electrons_per_configuration    = np.array([configuration_data.n_electrons   for configuration_data in self.configuration_data_list])

        # determine target_component for potential
        # For energy training with enhancement_factor_x: use raw v_x (symlog only, no LDA division)
        if include_potential:
            if target_component == "enhancement_factor_x":
                target_component_for_potential = "v_x"  # raw vxc, not v_x/v_lda
            else:
                target_component_for_potential = target_component.replace("e_", "v_")
        else:
            target_component_for_potential = None


        # Build configuration_data_list (one ExcConfiguration per configuration)
        X_list           : List[np.ndarray] = []
        y_energy_list    : List[np.ndarray] = []
        y_potential_list : List[np.ndarray] = []

        for configuration_data in self.configuration_data_list:
            # Extract features
            features_data_per_config = configuration_data.features_data.copy()

            # Extract target component for energy
            y_energy = configuration_data.extract_energy_target(
                target_component     = target_component,
                target_functional    = target_functional,
                reference_functional = reference_functional,
            )

            # Extract target component for potential, if needed
            if include_potential:
                y_potential = configuration_data.extract_potential_target(
                    target_component     = target_component_for_potential,
                    target_functional    = target_functional,
                    reference_functional = reference_functional,
                )
            else:
                y_potential = None

            # Store data
            X_list.append(features_data_per_config)
            y_energy_list.append(y_energy)
            y_potential_list.append(y_potential)

        # First symlog and reshape y once (shared by split and no-split paths)
        if use_symlog_for_features:
            X_list = [
                DataProcessor.symlog(features_data, linthresh=linthresh_for_features)
                for features_data in X_list
            ]
        if use_symlog_for_energy:
            y_energy_list = [
                DataProcessor.symlog(y_energy, linthresh=linthresh_for_energy)
                for y_energy in y_energy_list
            ]
        if include_potential and use_symlog_for_potential:
            y_potential_list = [
                DataProcessor.symlog(y_potential, linthresh=linthresh_for_potential)
                for y_potential in y_potential_list
            ]
        
        # reshape y if needed
        y_energy_list = [y_energy.reshape(-1, 1) for y_energy in y_energy_list]
        if include_potential:
            y_potential_list = [y_potential.reshape(-1, 1) for y_potential in y_potential_list]


        # Then compute scalers and weights data
        if split_train_val_test:
            
            # Split first, then fit scalers only on training set (by atomic number)
            unique_atoms = np.unique(atomic_numbers_per_configuration)
            
            # Ensure specified atoms are in training set
            if ensure_train_atoms is not None:
                ensure_train_atoms = np.array(ensure_train_atoms)
                atoms_guaranteed_train = np.intersect1d(unique_atoms, ensure_train_atoms)
                atoms_to_split = np.setdiff1d(unique_atoms, atoms_guaranteed_train)
                if len(atoms_guaranteed_train) > 0:
                    print(f"Ensuring atoms {atoms_guaranteed_train.tolist()} are in training set")
            else:
                atoms_guaranteed_train = np.array([], dtype=int)
                atoms_to_split = unique_atoms
            
            # Split remaining atoms into train/val/test
            if len(atoms_to_split) > 0:
                atoms_train_temp, atoms_temp = train_test_split(
                    atoms_to_split, test_size=(test_size + val_size), random_state=random_state
                )
                val_ratio = val_size / (test_size + val_size)
                atoms_val, atoms_test = train_test_split(
                    atoms_temp, test_size=(1 - val_ratio), random_state=random_state
                )
                # Combine guaranteed training atoms with randomly split training atoms
                atoms_train = np.concatenate([atoms_guaranteed_train, atoms_train_temp])
            else:
                # All atoms are guaranteed to be in training set
                atoms_train = atoms_guaranteed_train
                atoms_val   = np.array([], dtype=int)
                atoms_test  = np.array([], dtype=int)
            
            # Create masks (by atomic number) for train/val/test
            train_mask_per_configuration = np.isin(atomic_numbers_per_configuration, atoms_train)
            val_mask_per_configuration   = np.isin(atomic_numbers_per_configuration, atoms_val)
            test_mask_per_configuration  = np.isin(atomic_numbers_per_configuration, atoms_test)

            # Subsets for fit/transform (list indexing by boolean mask via zip)
            def _subset(lst, mask):
                return [x for x, m in zip(lst, mask) if m]
                
            X_train_list = _subset(X_list, train_mask_per_configuration)
            X_val_list   = _subset(X_list, val_mask_per_configuration)
            X_test_list  = _subset(X_list, test_mask_per_configuration)
            y_energy_train_list = _subset(y_energy_list, train_mask_per_configuration)
            y_energy_val_list   = _subset(y_energy_list, val_mask_per_configuration)
            y_energy_test_list  = _subset(y_energy_list, test_mask_per_configuration)
            if include_potential:
                y_potential_train_list = _subset(y_potential_list, train_mask_per_configuration)
                y_potential_val_list   = _subset(y_potential_list, val_mask_per_configuration)
                y_potential_test_list  = _subset(y_potential_list, test_mask_per_configuration)
            else:
                y_potential_train_list = None
                y_potential_val_list   = None
                y_potential_test_list  = None
            
            # Fit scalers on training set only; transform train/val/test (lists of ndarrays, possibly different lengths)
            scaler_X = None
            if scale_features:
                scaler_X = RobustScaler(**scaler_kwargs_for_features) if scaler_type_for_features == 'robust' else StandardScaler(**scaler_kwargs_for_features)
                scaler_X.fit(np.concatenate(X_train_list, axis=0))
                X_train_list = [scaler_X.transform(arr) for arr in X_train_list]
                X_val_list   = [scaler_X.transform(arr) for arr in X_val_list]
                X_test_list  = [scaler_X.transform(arr) for arr in X_test_list]

            scaler_y_energy = None
            if scale_energy:
                scaler_y_energy = RobustScaler(**scaler_kwargs_for_energy) if scaler_type_for_energy == 'robust' else StandardScaler(**scaler_kwargs_for_energy)
                scaler_y_energy.fit(np.concatenate(y_energy_train_list, axis=0))
                y_energy_train_list = [scaler_y_energy.transform(arr) for arr in y_energy_train_list]
                y_energy_val_list   = [scaler_y_energy.transform(arr) for arr in y_energy_val_list]
                y_energy_test_list  = [scaler_y_energy.transform(arr) for arr in y_energy_test_list]

            scaler_y_potential = None
            if include_potential and scale_potential:
                scaler_y_potential = RobustScaler(**scaler_kwargs_for_potential) if scaler_type_for_potential == 'robust' else StandardScaler(**scaler_kwargs_for_potential)
                scaler_y_potential.fit(np.concatenate(y_potential_train_list, axis=0))
                y_potential_train_list = [scaler_y_potential.transform(arr) for arr in y_potential_train_list]
                y_potential_val_list   = [scaler_y_potential.transform(arr) for arr in y_potential_val_list]
                y_potential_test_list  = [scaler_y_potential.transform(arr) for arr in y_potential_test_list]

        else:
            # Fit scalers on full data
            scaler_X = None
            if scale_features:
                scaler_X = RobustScaler(**scaler_kwargs_for_features) if scaler_type_for_features == 'robust' else StandardScaler(**scaler_kwargs_for_features)
                scaler_X.fit(np.concatenate(X_list, axis=0))
                X_list = [scaler_X.transform(arr) for arr in X_list]

            scaler_y_energy = None
            if scale_energy:
                scaler_y_energy = RobustScaler(**scaler_kwargs_for_energy) if scaler_type_for_energy == 'robust' else StandardScaler(**scaler_kwargs_for_energy)
                scaler_y_energy.fit(np.concatenate(y_energy_list, axis=0))
                y_energy_list = [scaler_y_energy.transform(arr) for arr in y_energy_list]
            
            scaler_y_potential = None
            if include_potential and scale_potential:
                scaler_y_potential = RobustScaler(**scaler_kwargs_for_potential) if scaler_type_for_potential == 'robust' else StandardScaler(**scaler_kwargs_for_potential)
                scaler_y_potential.fit(np.concatenate(y_potential_list, axis=0))
                y_potential_list = [scaler_y_potential.transform(arr) for arr in y_potential_list]
                    
        # Update weights data
        self.update_energy_weights_data(
            normalize_weight       = normalize_weight_for_energy,
            min_weight_ratio       = min_weight_ratio_for_energy,
            target_component       = target_component,
            target_functional      = target_functional,
            reference_functional   = reference_functional,
            linthresh_for_targets  = linthresh_for_energy,
            scaler_y               = scaler_y_energy,
            loss_norm              = loss_norm_for_energy,
            use_symlog_for_targets = use_symlog_for_energy,
        )

        if include_potential:
            self.update_potential_weights_data(
                normalize_weight       = normalize_weight_for_potential,
                min_weight_ratio       = min_weight_ratio_for_potential,
                target_component       = target_component_for_potential,
                target_functional      = target_functional,
                reference_functional   = reference_functional,
                linthresh_for_targets  = linthresh_for_potential,
                scaler_y               = scaler_y_potential,
                loss_norm              = loss_norm_for_potential,
                use_symlog_for_targets = use_symlog_for_potential,
            )
        

        # Build loader(s): one _make_loader that takes subset data + an/ne per config (use existing arrays, only expand to per-sample inside)
        def _make_loader(config_data_subset, X_sub, y_energy_sub, y_potential_sub, an_per_cfg, ne_per_cfg):
            configs: List[ExcConfiguration] = []

            for idx, configuration_data in enumerate(config_data_subset):

                # Compute total energy true
                E_true = configuration_data.compute_total_energy_true(
                    target_component     = target_component,
                    target_functional    = target_functional,
                    reference_functional = reference_functional,
                )
                # Obtain the number of elements and quadrature points
                dm = configuration_data.derivative_matrix if configuration_data.derivative_matrix is not None else self.shared_derivative_matrix
                n_elem = dm.shape[0]
                n_quad = dm.shape[1]

                # print out some debug information
                cfg = ExcConfiguration(
                    index                       = idx,
                    configuration_id            = configuration_data.configuration_id,
                    finite_element_number       = n_elem,
                    quadrature_point_number     = n_quad,
                    quadrature_nodes            = configuration_data.quadrature_nodes,
                    quadrature_weights          = configuration_data.quadrature_weights,
                    features                    = X_sub[idx],
                    features_physical           = configuration_data.features_data.copy(),
                    energy_density_transformed  = y_energy_sub[idx],
                    weights_for_energy_density  = configuration_data.energy_weights_data,
                    total_energy_true           = E_true,
                    potential_transformed       = (y_potential_sub[idx] if y_potential_sub else None) if include_potential else None,
                    weights_for_potential       = configuration_data.potential_weights_data if include_potential else None,
                    derivative_matrix           = configuration_data.derivative_matrix,
                )
                configs.append(cfg)

            return ExcDataLoader(
                configuration_data_list          = configs,
                include_potential                = include_potential,
                atomic_numbers_per_configuration = an_per_cfg,
                n_electrons_per_configuration    = ne_per_cfg,
                shared_derivative_matrix         = self.shared_derivative_matrix,
                target_functional                = target_functional,
                target_component                 = target_component,
                target_component_for_potential   = target_component_for_potential if include_potential else None,
                target_mode                      = target_mode,
                reference_functional             = reference_functional,
                features_list                    = self.features_list,
                scale_features                   = scale_features,
                scale_potential                  = scale_potential,
                scale_energy                     = scale_energy,
                scaler_type_for_features         = scaler_type_for_features,
                scaler_type_for_potential        = scaler_type_for_potential,
                scaler_type_for_energy           = scaler_type_for_energy,
                scaler_kwargs_for_features       = scaler_kwargs_for_features,
                scaler_kwargs_for_potential      = scaler_kwargs_for_potential,
                scaler_kwargs_for_energy         = scaler_kwargs_for_energy,
                scaler_for_features              = scaler_X,
                scaler_for_potential             = scaler_y_potential,
                scaler_for_energy                = scaler_y_energy,
                use_symlog_for_features          = use_symlog_for_features,
                use_symlog_for_potential         = use_symlog_for_potential,
                use_symlog_for_energy            = use_symlog_for_energy,
                linthresh_for_features           = linthresh_for_features,
                linthresh_for_potential          = linthresh_for_potential,
                linthresh_for_energy             = linthresh_for_energy,
                loss_norm_for_potential          = loss_norm_for_potential,
                loss_norm_for_energy             = loss_norm_for_energy,
            )

        if split_train_val_test:
            # Build loaders for train/val/test
            config_train = _subset(self.configuration_data_list, train_mask_per_configuration)
            config_val   = _subset(self.configuration_data_list, val_mask_per_configuration)
            config_test  = _subset(self.configuration_data_list, test_mask_per_configuration)
            an_train = atomic_numbers_per_configuration[train_mask_per_configuration]
            an_val   = atomic_numbers_per_configuration[val_mask_per_configuration]
            an_test  = atomic_numbers_per_configuration[test_mask_per_configuration]
            ne_train = n_electrons_per_configuration[train_mask_per_configuration]
            ne_val   = n_electrons_per_configuration[val_mask_per_configuration]
            ne_test  = n_electrons_per_configuration[test_mask_per_configuration]
            return (
                _make_loader(config_train, X_train_list, y_energy_train_list, y_potential_train_list if include_potential else None, an_train, ne_train),
                _make_loader(config_val,   X_val_list,   y_energy_val_list,   y_potential_val_list   if include_potential else None, an_val,   ne_val),
                _make_loader(config_test,  X_test_list,  y_energy_test_list,  y_potential_test_list  if include_potential else None, an_test,  ne_test),
            )

        # Build loader for full data
        return _make_loader(
            config_data_subset = self.configuration_data_list, 
            X_sub              = X_list, 
            y_energy_sub       = y_energy_list, 
            y_potential_sub    = y_potential_list if include_potential else None, 
            an_per_cfg         = atomic_numbers_per_configuration, 
            ne_per_cfg         = n_electrons_per_configuration
        )


    def update_potential_weights_data(
        self, 
        normalize_weight      : bool = True,  # If True, use N_v^(c) normalization; if False, N_v^(c) = 1
        min_weight_ratio      : float = 1e-3,
        target_component      : Optional[PotentialTargetComponent] = None,
        target_functional     : Optional[str] = None,
        reference_functional  : Optional[str] = None,
        linthresh_for_targets : Optional[float] = None,
        scaler_y              : ScalerType = None,
        loss_norm             : LossNorm = "L1norm",
        use_symlog_for_targets: bool = True,
    ):
        """
        Update the potential weights for weighted per-grid losses on transformed targets.

        When ``use_symlog_for_targets`` is True (potential targets are symlogged before scaling),
        Jacobian factors use ``|Symlog'[Δv_xc]|`` as in the standard formulas.

        When False (no symlog on targets, only optional robust/standard scaling), treat
        ``|Symlog'| ≡ 1`` so weights are ρ–quadrature-consistent with scaled raw Δv only.

        **L1norm** (MAE-style; default): ``w_i ∝ (4π r_i² σ_v w ρ N_v^(c)) / |Symlog'|`` if symlog,
        else ``w_i ∝ 4π r_i² σ_v w ρ N_v^(c)``.

        **L2norm** (MSE-style): extra factor ``σ_v / |Symlog'|`` with symlog; without symlog,
        extra factor ``σ_v`` only (equivalent to ``|Symlog'|=1``).

        If normalize_weight=True, N_v^(c) normalizes each configuration so all atoms contribute equally
        to the total weight sum; if False, N_v^(c)=1.
        """
        _validate_loss_norm("loss_norm", loss_norm)
        # Check if necessary parameters are provided
        assert target_component is not None, \
            TARGET_COMPONENT_NOT_PROVIDED_ERROR.format("symlog")
        assert target_functional is not None, \
            TARGET_FUNCTIONAL_NOT_PROVIDED_ERROR.format("symlog")
        if use_symlog_for_targets:
            assert linthresh_for_targets is not None, \
                LINTHRESH_FOR_TARGETS_NOT_PROVIDED_ERROR.format("symlog")

        # check if weights data is already updated
        if self.potential_weights_data_is_updated:
            print(POTENTIAL_WEIGHTS_DATA_ALREADY_UPDATED_WARNING)

        # Get σ_v from scaler_y
        sigma_v = 1.0

        if scaler_y is not None:
            if hasattr(scaler_y, 'scale_'):
                # StandardScaler or RobustScaler
                sigma_v = scaler_y.scale_
                # Convert to scalar if it's an array
                if isinstance(sigma_v, np.ndarray):
                    if sigma_v.ndim > 0:
                        sigma_v = float(sigma_v[0])  # Take first component if multi-dimensional
                    else:
                        sigma_v = float(sigma_v)
                else:
                    sigma_v = float(sigma_v)
        

        # Compute weights data for each configuration
        raw_weights_data_list = []
        normalization_sum_list = []  # For computing N_v^(c) if needed

        for configuration_data in self.configuration_data_list:
            rho = configuration_data.rho_filtered.reshape(-1, 1)
            quadrature_nodes = configuration_data.quadrature_nodes_filtered.reshape(-1, 1)
            
            # Get quadrature weights for this configuration
            # Default behavior: use quadrature_weights_filtered from configuration_data
            if hasattr(configuration_data, "quadrature_weights_filtered"):
                w = configuration_data.quadrature_weights_filtered.reshape(-1, 1)
            else:
                # Fallback: uniform weights if quadrature weights are not available
                w = np.ones_like(quadrature_nodes)
            
            # Compute Δv_xc^(c)(r)
            delta_v_xc = configuration_data.extract_potential_target(
                target_component     = target_component,
                target_functional    = target_functional,
                reference_functional = reference_functional,
            )
            
            # Ensure 2D shape
            if delta_v_xc.ndim == 1:
                delta_v_xc = delta_v_xc.reshape(-1, 1)
            if rho.ndim == 1:
                rho = rho.reshape(-1, 1)

            if use_symlog_for_targets:
                abs_symlog_derivative = np.asarray(
                    DataProcessor.symlog_derivative(
                        delta_v_xc, linthresh=linthresh_for_targets
                    ),
                    dtype=float,
                )
            else:
                abs_symlog_derivative = np.ones_like(delta_v_xc, dtype=float)
            safe_dsym = np.maximum(abs_symlog_derivative, 1e-300)

            # Calculate numerator: 4πr²σ_v w(r)ρ(r)
            numerator = 4 * np.pi * quadrature_nodes**2 * sigma_v * w * rho
            if target_component == "enhancement_factor_x":
                from ..xc.lda import lda_exchange_potential_generic
                lda_exchange_potential = lda_exchange_potential_generic(rho)
                numerator = numerator * np.abs(lda_exchange_potential)

            l1_weights = numerator / abs_symlog_derivative
            if loss_norm == "L2norm":
                weights_data = l1_weights * (sigma_v / safe_dsym)
            else:
                weights_data = l1_weights

            # Ensure non-negative (allow very small values; actual lower bound set globally later)
            weights_data = np.maximum(weights_data, 0.0)

            raw_weights_data_list.append(weights_data)

            # If using normalization, compute the sum for N_v^(c)
            if normalize_weight:
                normalization_sum = float(np.sum(weights_data))
                normalization_sum_list.append(normalization_sum)
            else:
                normalization_sum_list.append(None)

        # Apply normalization factor N_v^(c) if needed
        if normalize_weight:
            normalized_weights_data_list = []
            for weights_data, norm_sum in zip(raw_weights_data_list, normalization_sum_list):
                if norm_sum is not None and norm_sum > 0:
                    N_v_c = 1.0 / norm_sum
                    normalized_weights = weights_data * N_v_c
                else:
                    normalized_weights = weights_data
                normalized_weights_data_list.append(normalized_weights)
        else:
            normalized_weights_data_list = raw_weights_data_list

        # Apply minimum weight ratio constraint (global lower bound across all configurations)
        final_weights_data_list = []
        if len(normalized_weights_data_list) > 0:
            # Compute global maximum weight over all configurations
            global_max_weight = max(
                float(np.max(normalized_weights_data))
                for normalized_weights_data in normalized_weights_data_list
            )
            # Global lower bound controlled by min_weight_ratio
            tiny = 1e-12  # purely for numerical safety
            global_min_weight = max(min_weight_ratio * global_max_weight, tiny)
            
            # Debug: Print global bounds info
            if len(normalized_weights_data_list) > 1:
                print(f"Global potential weights bound: min={global_min_weight:.6e}, max={global_max_weight:.6e}, ratio={min_weight_ratio:.6e}")
        else:
            global_min_weight = 0.0

        for normalized_weights_data in normalized_weights_data_list:
            # Apply global minimum weight constraint
            final_weights = np.maximum(
                normalized_weights_data,
                global_min_weight
            )
            final_weights_data_list.append(final_weights)

        # set the weights data
        for configuration_data, final_weights_data in \
            zip(self.configuration_data_list, final_weights_data_list):
            configuration_data.set_potential_weights_data(final_weights_data)

        # update the weights data is updated flag
        self.potential_weights_data_is_updated = True



    def update_energy_weights_data(
        self, 
        normalize_weight      : bool = True,  # If True, use N_e^(c) normalization; if False, N_e^(c) = 1
        min_weight_ratio      : float = 1e-3,
        target_component      : Optional[EnergyTargetComponent] = None,
        target_functional     : Optional[str] = None,
        reference_functional  : Optional[str] = None,
        linthresh_for_targets : Optional[float] = None,
        scaler_y              : ScalerType = None,
        loss_norm             : LossNorm = "L1norm",
        use_symlog_for_targets: bool = True,
    ):
        """
        Per-grid weights for energy density training (analogous to potential).

        With ``use_symlog_for_targets=True``, weights include ``1/|Symlog'[Δe_xc]|`` (and L2norm uses
        ``σ_e/|Symlog'|`` as for potential). With ``False``, use ``|Symlog'|≡1`` (scaled physical
        targets only).
        """
        _validate_loss_norm("loss_norm", loss_norm)
        # Check if necessary parameters are provided
        assert target_component is not None, \
            TARGET_COMPONENT_NOT_PROVIDED_ERROR.format("symlog")
        assert target_functional is not None, \
            TARGET_FUNCTIONAL_NOT_PROVIDED_ERROR.format("symlog")
        if use_symlog_for_targets:
            assert linthresh_for_targets is not None, \
                LINTHRESH_FOR_TARGETS_NOT_PROVIDED_ERROR.format("symlog")

        # check if weights data is already updated
        if self.energy_weights_data_is_updated:
            print(ENERGY_WEIGHTS_DATA_ALREADY_UPDATED_WARNING)

        # Get σ_e from scaler_y
        sigma_e = 1.0

        if scaler_y is not None:
            if hasattr(scaler_y, 'scale_'):
                # StandardScaler or RobustScaler
                sigma_e = scaler_y.scale_
                # Convert to scalar if it's an array
                if isinstance(sigma_e, np.ndarray):
                    if sigma_e.ndim > 0:
                        sigma_e = float(sigma_e[0])  # Take first component if multi-dimensional
                    else:
                        sigma_e = float(sigma_e)
                else:
                    sigma_e = float(sigma_e)

        # Compute weights data for each configuration
        raw_weights_data_list = []
        normalization_sum_list = []  # For computing N_e^(c) if needed

        for configuration_data in self.configuration_data_list:
            rho = configuration_data.rho_filtered.reshape(-1, 1)
            quadrature_nodes = configuration_data.quadrature_nodes_filtered.reshape(-1, 1)

            # Get quadrature weights for this configuration
            if hasattr(configuration_data, "quadrature_weights_filtered"):
                w = configuration_data.quadrature_weights_filtered.reshape(-1, 1)
            else:
                w = np.ones_like(quadrature_nodes)

            # Compute Δe_xc^(c)(r)
            delta_e_xc = configuration_data.extract_energy_target(
                target_component     = target_component,
                target_functional    = target_functional,
                reference_functional = reference_functional,
            )

            # Ensure 2D shape
            if delta_e_xc.ndim == 1:
                delta_e_xc = delta_e_xc.reshape(-1, 1)

            if use_symlog_for_targets:
                abs_symlog_derivative = np.asarray(
                    DataProcessor.symlog_derivative(
                        delta_e_xc, linthresh=linthresh_for_targets
                    ),
                    dtype=float,
                )
            else:
                abs_symlog_derivative = np.ones_like(delta_e_xc, dtype=float)
            safe_dsym = np.maximum(abs_symlog_derivative, 1e-300)

            # Calculate numerator: 4πr²σ_e w(r)ρ(r)
            numerator = 4 * np.pi * quadrature_nodes**2 * sigma_e * w * rho
            if target_component == "enhancement_factor_x":
                from ..xc.lda import lda_exchange_energy_density_generic
                lda_exchange_energy_density = lda_exchange_energy_density_generic(rho) * rho
                numerator = numerator * np.abs(lda_exchange_energy_density)

            l1_weights = numerator / abs_symlog_derivative
            if loss_norm == "L2norm":
                weights_data = l1_weights * (sigma_e / safe_dsym)
            else:
                weights_data = l1_weights

            # Ensure non-negative (allow very small values; actual lower bound set globally later)
            weights_data = np.maximum(weights_data, 0.0)

            raw_weights_data_list.append(weights_data)

            # If using normalization, compute the sum for N_e^(c)
            if normalize_weight:
                normalization_sum = float(np.sum(weights_data))
                normalization_sum_list.append(normalization_sum)
            else:
                normalization_sum_list.append(None)

        # Apply normalization factor N_e^(c) if needed
        if normalize_weight:
            normalized_weights_data_list = []
            for weights_data, norm_sum in zip(raw_weights_data_list, normalization_sum_list):
                if norm_sum is not None and norm_sum > 0:
                    N_e_c = 1.0 / norm_sum
                    normalized_weights = weights_data * N_e_c
                else:
                    normalized_weights = weights_data
                normalized_weights_data_list.append(normalized_weights)
        else:
            normalized_weights_data_list = raw_weights_data_list

        # Apply minimum weight ratio constraint (global lower bound across all configurations)
        final_weights_data_list = []
        if len(normalized_weights_data_list) > 0:
            global_max_weight = max(
                float(np.max(normalized_weights_data))
                for normalized_weights_data in normalized_weights_data_list
            )
            tiny = 1e-12  # purely for numerical safety
            global_min_weight = max(min_weight_ratio * global_max_weight, tiny)

            if len(normalized_weights_data_list) > 1:
                print(f"Global energy density weights bound: min={global_min_weight:.6e}, max={global_max_weight:.6e}, ratio={min_weight_ratio:.6e}")
        else:
            global_min_weight = 0.0

        for normalized_weights_data in normalized_weights_data_list:
            final_weights = np.maximum(
                normalized_weights_data,
                global_min_weight
            )
            final_weights_data_list.append(final_weights)

        # set the weights data
        for configuration_data, final_weights_data in \
            zip(self.configuration_data_list, final_weights_data_list):
            configuration_data.set_energy_weights_data(final_weights_data)

        # update the weights data is updated flag
        self.energy_weights_data_is_updated = True


    @staticmethod
    def extract_potential_component(xc_data: XCDataType, target_component: PotentialTargetComponent) -> np.ndarray:
        """
        Extract the component of the XC data (potential: v_x, v_c, v_xc, v_x_v_c).
        """
        return SingleConfigurationData.extract_potential_component(xc_data, target_component)


    @staticmethod
    def extract_energy_component(xc_data: XCDataType, target_component: EnergyTargetComponent) -> np.ndarray:
        """
        Extract the energy density component from XC data.
        target_component: "e_xc", "e_x", "e_c", or "e_x_e_c".
        Returns 1d array (e_xc/e_x/e_c) or 2d array (e_x_e_c).
        """
        return SingleConfigurationData.extract_energy_component(xc_data, target_component)


    @staticmethod
    def integrate_energy_density(
        quadrature_nodes   : np.ndarray,
        quadrature_weights : np.ndarray,
        energy_density     : np.ndarray,
    ) -> np.ndarray:
        """
        Integrate energy density over radial grid to get total energy.
        E = ∫ e(r) * 4πr² dr ≈ ∑_i (4π * r_i² * w_i * e(r_i)).
        energy_density: shape (n_points,) or (n_points, n_targets). Returns scalar or (n_targets,).
        """
        r = np.asarray(quadrature_nodes).reshape(-1, 1)
        w = np.asarray(quadrature_weights).reshape(-1, 1)
        e = np.asarray(energy_density)
        if e.ndim == 1:
            e = e.reshape(-1, 1)
        integrand = 4 * np.pi * (r ** 2) * w * e  # (n_points, n_targets)
        return np.sum(integrand, axis=0)


    @staticmethod
    def calculate_symlog_weights(
        rho                 : np.ndarray,
        vxc_target_physical : np.ndarray,
        linthresh           : float,
        min_weight          : float = 1e-6
    ) -> np.ndarray:
        """
        Calculate symlog-based weights from density and target Vxc values.
        """
        return DataProcessor.calculate_symlog_weights(rho, vxc_target_physical, linthresh, min_weight)


    @staticmethod
    def normalize_weights_by_atom(
        weights          : np.ndarray,
        atomic_numbers   : np.ndarray,
        min_weight_ratio : float = 1e-2
    ) -> np.ndarray:
        """
        Normalize weights per atom to keep relative scales comparable.
        """
        raise NotImplementedError("This method is deprecated. Use method self.update_potential_weights_data() instead.")
        
        return DataProcessor.normalize_weights_by_atom(weights, atomic_numbers, min_weight_ratio)


    @staticmethod
    def normalize_weights_by_configuration(
        raw_weights_data_list         : List[np.ndarray],
        weights_sum_per_configuration : float = 100.0,
        min_weight_ratio              : float = 1e-2,
    ) -> List[np.ndarray]:
        """
        Normalize weights per configuration to keep relative scales comparable.
        Here, we let the weights for each configuration to have the same sum.
        """
        normalized_weights_data_list = []

        for raw_weights_data in raw_weights_data_list:
            # Let the weights for each configuration to have the same sum.
            normalized_weights_data = raw_weights_data / np.sum(raw_weights_data) * weights_sum_per_configuration

            # apply minimum weight constraint
            normalized_weights_data = np.maximum(normalized_weights_data, min_weight_ratio * np.max(normalized_weights_data))
            normalized_weights_data_list.append(normalized_weights_data)

        return normalized_weights_data_list


    def get_scf_folder_path(self, atomic_number: int) -> str:
        """
        Get SCF folder path for a given atomic number.
        Note: This method finds the configuration folder by matching atomic_number in meta.json files.
        It supports both new format (configuration_XXX) and old format (atom_XXX) for backward compatibility.
        """
        assert isinstance(atomic_number, int), \
            ATOMIC_NUMBER_NOT_INTEGER_ERROR.format(atomic_number)
        assert atomic_number in self.atomic_numbers_unique, \
            ATOMIC_NUMBER_NOT_IN_DATASET_ERROR.format(atomic_number)
        
        # Find configuration folder by searching for matching atomic_number in meta.json
        # This handles both new and old folder naming formats
        for config_data in self.configuration_data_list:
            if config_data.atomic_number == atomic_number:
                # Use the folder path from the loaded configuration data
                return os.path.join(config_data.scf_folder_path)
        
        # Fallback: try to find by searching all folders (for backward compatibility)
        # This should rarely be needed if configuration_data_list is properly loaded
        import json
        for item in os.listdir(self.data_root):
            item_path = os.path.join(self.data_root, item)
            if not os.path.isdir(item_path):
                continue
            if item.startswith('configuration_') or item.startswith('atom_'):
                meta_path = os.path.join(item_path, "meta.json")
                if os.path.exists(meta_path):
                    try:
                        with open(meta_path, "r") as meta_file:
                            meta_data = json.load(meta_file)
                        if int(meta_data.get("atomic_number", 0)) == atomic_number:
                            return os.path.join(item_path, self.scf_xc_functional.lower())
                    except Exception:
                        continue
        
        # If still not found, raise error
        raise ValueError(ATOMIC_NUMBER_NOT_IN_DATASET_ERROR.format(atomic_number))


    def exists_functional(self, functional: str) -> bool:
        """
        Check if a functional exists in the dataset.
        """
        return functional == self.scf_xc_functional or functional in self.forward_pass_xc_functional_list


    def get_xc_data(self, functional: str) -> XCDataType:
        """
        Get XC data for a given functional.
        """
        if functional == self.scf_xc_functional:
            return self.scf_xc_data
        if functional in self.forward_pass_xc_functional_list:
            return self.forward_pass_xc_data_list[self.forward_pass_xc_functional_list.index(functional)]
        raise ValueError(FUNCTIONAL_NOT_IN_DATASET_ERROR.format(functional, self.scf_xc_functional, self.forward_pass_xc_functional_list))


    def get_features_data(self, feature_name: str) -> np.ndarray:
        """
        Get features data for a given feature name.
        """
        assert isinstance(feature_name, str),\
            FEATURE_NAME_NOT_STRING_ERROR.format(feature_name)
        assert feature_name in self.features_list,\
            FEATURE_NAME_NOT_IN_FEATURES_LIST_ERROR.format(feature_name, self.features_list)
        return self.features_data[:, self.features_list.index(feature_name)]


    @property
    def rho(self) -> np.ndarray:
        """
        Get density data for all configurations.
        """
        assert "rho" in self.features_list,\
            FEATURE_NAME_NOT_IN_FEATURES_LIST_ERROR.format("rho", self.features_list)
        return self.features_data[:, self.features_list.index("rho")]


    @property
    def potential_weights_data(self) -> np.ndarray:
        """
        Get weights data for all configurations.
        """
        assert self.potential_weights_data_is_updated
        return np.concatenate([configuration_data.potential_weights_data for configuration_data in self.configuration_data_list], axis=0)


    @property
    def energy_weights_data(self) -> np.ndarray:
        """
        Get energy weights data for all configurations.
        """
        assert self.energy_weights_data_is_updated
        return np.concatenate([configuration_data.energy_weights_data for configuration_data in self.configuration_data_list], axis=0)


    @property
    def atomic_numbers_unique(self) -> np.ndarray:
        """
        Get unique atomic numbers.
        """
        if self.cached_atomic_numbers_unique is None:
            self.cached_atomic_numbers_unique = np.sort(np.unique(self.atomic_numbers_per_configuration))
        return self.cached_atomic_numbers_unique
    

    @property
    def atomic_numbers_per_configuration(self) -> np.ndarray:
        """
        Get atomic numbers per configuration.
        """
        if self.cached_atomic_numbers_per_configuration is None:
            self.cached_atomic_numbers_per_configuration = np.array([configuration_data.atomic_number for configuration_data in self.configuration_data_list])
        return self.cached_atomic_numbers_per_configuration
    

    @property
    def atomic_numbers_per_sample(self) -> np.ndarray:
        """
        Get atomic numbers per sample.
        """
        if self.cached_atomic_numbers_per_sample is None:
            self.cached_atomic_numbers_per_sample = np.concatenate([configuration_data.atomic_numbers for configuration_data in self.configuration_data_list])
        return self.cached_atomic_numbers_per_sample


    @property
    def n_electrons_per_sample(self) -> np.ndarray:
        """
        Get number of electrons per sample.
        """
        if self.cached_n_electrons_per_sample is None:
            self.cached_n_electrons_per_sample = np.concatenate([np.full(configuration_data.n_filtered, configuration_data.n_electrons) for configuration_data in self.configuration_data_list])
        return self.cached_n_electrons_per_sample


    @property
    def n_atoms(self) -> int:
        """
        Get number of atomic numbers.
        """
        return len(self.atomic_numbers_unique)
    

    @property
    def n_configurations(self) -> int:
        """
        Get number of configurations.
        """
        return len(self.configuration_data_list)
    

    @property
    def n_samples(self) -> int:
        """
        Get number of samples.
        """
        return len(self.atomic_numbers_per_sample)
    

    @property
    def cutoff_radii(self) -> np.ndarray:
        """
        Get cutoff radii.
        """
        if self.cached_cutoff_radii is None:
            self.cached_cutoff_radii = np.array([configuration_data.cutoff_radius for configuration_data in self.configuration_data_list])
        return self.cached_cutoff_radii
    

    @property
    def cutoff_indices(self) -> np.ndarray:
        """
        Get cutoff indices.
        """
        if self.cached_cutoff_indices is None:
            self.cached_cutoff_indices = np.array([configuration_data.cutoff_idx for configuration_data in self.configuration_data_list])
        return self.cached_cutoff_indices


    @property
    def quadrature_nodes(self) -> np.ndarray:
        """
        Get quadrature nodes for all configurations.
        """
        if self.cached_quadrature_nodes is None:
            self.cached_quadrature_nodes = np.concatenate([configuration_data.quadrature_nodes_filtered for configuration_data in self.configuration_data_list])
        return self.cached_quadrature_nodes


    @property
    def configuration_ids(self) -> List[int]:
        """
        Get configuration ids for all configurations.
        """
        if self.cached_configuration_ids is None:
            self.cached_configuration_ids = [configuration_data.configuration_id for configuration_data in self.configuration_data_list]
        return self.cached_configuration_ids


    @property
    def features_data(self) -> np.ndarray:
        """
        Get features data for all configurations.
        """
        if self.cached_features_data is None:
            self.cached_features_data = np.concatenate([configuration_data.features_data for configuration_data in self.configuration_data_list])
        return self.cached_features_data


    @property
    def scf_xc_data(self) -> XCDataType:
        """
        Get SCF XC data for all configurations.
        """
        if self.cached_scf_xc_data is None:
            scf_xc_data_list = [configuration_data.scf_xc_data for configuration_data in self.configuration_data_list]
            self.cached_scf_xc_data = (
                np.concatenate([v_x for v_x, _, _, _ in scf_xc_data_list]),
                np.concatenate([v_c for _, v_c, _, _ in scf_xc_data_list]),
                np.concatenate([e_x for _, _, e_x, _ in scf_xc_data_list]) if self.include_energy_density else None,
                np.concatenate([e_c for _, _, _, e_c in scf_xc_data_list]) if self.include_energy_density else None,
            )
        return self.cached_scf_xc_data
    

    @property
    def forward_pass_xc_data_list(self) -> List[XCDataType]:
        """
        Get forward pass XC data for all configurations.
        """
        if self.cached_forward_pass_xc_data_list is None:

            # nest forward pass XC data into lists
            forward_pass_xc_data_list_list = [[] for _ in range(len(self.forward_pass_xc_functional_list))]
            for configuration_data in self.configuration_data_list:
                for idx, forward_pass_xc_data in enumerate(configuration_data.forward_pass_xc_data_list):
                    forward_pass_xc_data_list_list[idx].append(forward_pass_xc_data)

            # convert lists to numpy arrays
            self.cached_forward_pass_xc_data_list = [
                (
                    np.concatenate([v_x for v_x, _, _, _ in forward_pass_xc_data_list]),
                    np.concatenate([v_c for _, v_c, _, _ in forward_pass_xc_data_list]),
                    np.concatenate([e_x for _, _, e_x, _ in forward_pass_xc_data_list]) if self.include_energy_density else None,
                    np.concatenate([e_c for _, _, _, e_c in forward_pass_xc_data_list]) if self.include_energy_density else None,
                )
                for forward_pass_xc_data_list in forward_pass_xc_data_list_list
            ]
        return self.cached_forward_pass_xc_data_list





class AtomicDataManager:
    
    """
    Main class for managing atomic data, consisting of:
    - DataGenerator : generating atomic data
    - DataLoader    : loading atomic data
    - DataProcessor : processing atomic data
    """

    def __init__(self,
        data_root                   : str,
        scf_xc_functional           : str,
        forward_pass_xc_functionals : Optional[str | list[str]],
        auto_confirm                : bool = False,  # If True, automatically confirm prompts without user input
    ):

        """
        Args:
            data_root                   : Root directory of the dataset
            scf_xc_functional           : SCF XC functional (used for full SCF calculation to convergence)
            forward_pass_xc_functionals : Forward pass XC functional(s), if not None, will perform forward pass for each functional based on SCF results
            auto_confirm                : If True, automatically confirm prompts without user input (useful for batch/background jobs). Defaults to False.
        """

        # Check if data root exists
        if not os.path.exists(data_root):
            # Safety confirmation prompt for creating empty directory
            print("\n" + "=" * 60)
            print("WARNING: data root '{}' does not found.".format(data_root))
            print("=" * 60)
            if auto_confirm:
                print("\nAuto-confirming: Creating directory...")
                os.makedirs(data_root)
            else:
                if input("\nDo you want to create an empty directory? (y/n): ").strip().lower() == 'y':
                    os.makedirs(data_root)
                else:
                    print("Operation cancelled by user.")
                    exit(0)
            
            print("Directory created successfully.")


        # Process forward_pass_xc_functionals: None -> [], str -> [str], list -> list
        if forward_pass_xc_functionals is None:
            forward_pass_xc_functional_list = []
        elif isinstance(forward_pass_xc_functionals, str):
            forward_pass_xc_functional_list = [forward_pass_xc_functionals]
        elif isinstance(forward_pass_xc_functionals, list):
            forward_pass_xc_functional_list = forward_pass_xc_functionals
        else:
            raise TypeError(FORWARD_PASS_XC_FUNCTIONAL_NOT_NONE_STR_OR_LIST_ERROR.format(type(forward_pass_xc_functionals)))


        # Type and value checks
        self.check_data_root(data_root)
        self.check_scf_xc_functional(scf_xc_functional)
        self.check_forward_pass_xc_functional_list(forward_pass_xc_functional_list)


        # Initialize attributes
        self.data_root                       : str       = data_root
        self.scf_xc_functional               : str       = scf_xc_functional
        self.forward_pass_xc_functional_list : List[str] = forward_pass_xc_functional_list
        self.auto_confirm                    : bool      = auto_confirm

        # Initialize tools
        self.generator : DataGenerator = DataGenerator()
        self.loader    : DataLoader    = DataLoader()
        self.processor : DataProcessor = DataProcessor()

        # Initialize default parameters
        ## For data generation
        self.default_radius_cutoff_rho_threshold : float = 1e-6
        self.default_radius_cutoff_v_x_threshold : float = 1e-8
        self.default_radius_cutoff_v_c_threshold : float = 1e-8

        ## For data processing
        self.default_smooth_radius_threshold : float = 5.0
        self.default_smooth_method           : str = 'savgol'
        self.default_smooth_kwargs           : Dict = {}

        ## Log-linear asymptotic rho smoothing (``DataProcessor.smooth_rho_with_loglinear_asymptote``, used in ``load_data``)
        self.default_smooth_rho_r_a       : float = 7.0
        self.default_smooth_rho_r_b       : float = 10.0
        self.default_smooth_rho_r_c       : float = 13.0
        self.default_smooth_rho_log_floor : float = 1e-30


    # Data generation methods
    def generate_data(self, 
        # Required arguments
        atomic_number_list        : list[int | float], 
        n_electrons_list          : Optional[list[int | float]] = None,
        use_oep                   : Optional[bool] = None,  # If None, uses default from XC_FUNCTIONAL_OEP_DEFAULT based on scf_xc_functional

        # Optional arguments controlling the contents of the dataset
        save_energy_density       : bool  = False,
        save_intermediate         : bool  = False,
        save_full_spectrum        : bool  = False,
        save_full_orbitals        : bool  = False,
        save_derivative_matrix    : bool  = False,
        start_configuration_index : int   = 1,

        # Optional arguments controlling the generation process
        # Grid, basis, and mesh parameters
        domain_size               : float = 20.0,
        finite_elements_number    : int   = 35,
        polynomial_order          : int   = 20,
        quadrature_point_number   : int   = 43,
        oep_basis_number          : int   = 5,
        mesh_type                 : str   = "polynomial",
        mesh_concentration        : float = 2.0,
        mesh_spacing              : float = 0.1,

        # SCF convergence parameters
        scf_tolerance             : float = 1e-8,
        max_scf_iterations        : int   = 500,
        max_scf_iterations_outer  : Optional[int] = None,
        use_pulay_mixing          : bool  = True,
        use_preconditioner        : bool  = True,
        pulay_mixing_parameter    : float = 1.0,
        pulay_mixing_history      : int   = 7,
        pulay_mixing_frequency    : int   = 3,
        linear_mixing_alpha1      : float = 0.75,
        linear_mixing_alpha2      : float = 0.95,

        # Advanced functional parameters
        hybrid_mixing_parameter           : float = None,
        frequency_quadrature_point_number : int   = None,
        angular_momentum_cutoff           : int   = None,
        double_hybrid_flag                : bool  = None,
        oep_mixing_parameter              : float = None,
        enable_parallelization            : bool  = None,

        # Debugging and verbose parameters
        verbose                   : bool  = True,
        overwrite                 : Optional[bool] = None,  # If True, automatically confirm prompts without user input. If None, uses instance default.

        # deprecated arguments
        finite_elements           : Optional[int] = None,
    ):
        """
        Generate atomic dataset, based on SCF and forward pass XC functionals.
        

        Required arguments
        ------------------
        `atomic_number_list` : list[int | float]
            List of atomic numbers to generate data for (can be fractional).
        `n_electrons_list` : list[int | float] | None
            List of number of electrons to generate data for (can be fractional). If None, defaults to atomic_number_list.

        Dataset content control
        -----------------------
        `save_energy_density` : bool
            Whether to save energy density in the dataset. Defaults to False.
        `save_intermediate` : bool
            Whether to save intermediate information during SCF. Defaults to False.
        `save_full_spectrum` : bool
            Whether to save full_eigen_energies and full_l_terms. Defaults to False.
        `save_full_orbitals` : bool
            Whether to save full_orbitals (default: False, to save storage).
        `save_derivative_matrix` : bool
            Whether to save derivative matrix. Most systems have the same derivative matrix when using
            the same grid/basis/mesh parameters, so a shared derivative matrix is saved at the dataset root.
            If an atom's derivative matrix differs from the shared one, it is saved locally and recorded in meta.json.
            Defaults to False.
        `start_configuration_index` : int
            Starting configuration index for generated folders. Configuration folders will be named
            configuration_XXX starting from this index. Default is 1.
            For example, if start_configuration_index=5, the first atom will be saved as configuration_005,
            the second as configuration_006, etc.

        Grid, basis, and mesh parameters
        --------------------------------
        `domain_size` : float
            Radial computational domain size in atomic units (typically 10-30 Bohr). Defaults to 20.0.
        `finite_elements_number` : int
            Number of finite elements in the computational domain. Defaults to 35.
        `polynomial_order` : int
            Polynomial order of basis functions within each finite element (typically 20-40). Defaults to 20.
        `quadrature_point_number` : int
            Number of quadrature points for numerical integration (recommended: 3-4x polynomial_order). Defaults to 43.
        `oep_basis_number` : int
            Basis size used in OEP calculations when enabled. Defaults to 5.
        `mesh_type` : str
            Mesh distribution type ('exponential', 'polynomial', 'uniform'). Defaults to 'polynomial'.
        `mesh_concentration` : float
            Mesh concentration parameter (controls point density distribution). Defaults to 2.0.
        `mesh_spacing` : float
            Used to set the output uniform mesh spacing, irrelevant during SCF calculation. Defaults to 0.1.

        Self-consistent field (SCF) convergence parameters
        --------------------------------------------------
        `scf_tolerance` : float
            SCF convergence tolerance (typically 1e-8). Defaults to 1e-8 (1e-6 for SCAN/RSCAN/R2SCAN functionals).
        `max_scf_iterations` : int
            Maximum number of inner SCF iterations. If None, uses default (500). Defaults to None.
        `max_scf_iterations_outer` : int | None
            Maximum number of outer SCF iterations (for functionals requiring outer loop like HF, EXX, RPA, PBE0).
            If None, uses default (50 when needed, otherwise not used). Defaults to None.
        `use_pulay_mixing` : bool
            True for Pulay mixing for SCF convergence, False for linear mixing. Defaults to True.
        `use_preconditioner` : bool
            Flag for using preconditioner for SCF convergence. Defaults to True.
        `pulay_mixing_parameter` : float
            Pulay mixing parameter. Defaults to 1.0.
        `pulay_mixing_history` : int
            Pulay mixing history. Defaults to 7.
        `pulay_mixing_frequency` : int
            Pulay mixing frequency. Defaults to 3.
        `linear_mixing_alpha1` : float
            Linear mixing parameter (alpha_1). Defaults to 0.75.
        `linear_mixing_alpha2` : float
            Linear mixing parameter (alpha_2). Defaults to 0.95.

        Advanced functional parameters
        ------------------------------
        `hybrid_mixing_parameter` : float
            Mixing parameter for hybrid/double-hybrid functionals. Defaults to 1.0.
        `frequency_quadrature_point_number` : int
            Number of frequency quadrature points for RPA calculations. Defaults to 25.
        `angular_momentum_cutoff` : int
            Maximum angular momentum quantum number to include. Defaults to 4.
        `double_hybrid_flag` : bool
            Flag for double-hybrid functional methods. Defaults to False.
        `oep_mixing_parameter` : float
            Scaling parameter (λ) for OEP exchange/correlation potentials. Defaults to 1.0.
        `enable_parallelization` : bool
            Flag for parallelization of RPA calculations. Defaults to False.

        Debugging and verbose parameters
        --------------------------------
        `verbose` : bool
            Whether to print information during execution. Defaults to True.
        `overwrite` : bool | None
            If True, automatically confirm prompts without user input (useful for batch/background jobs).
            If None, uses the value set in AtomicDataManager.__init__(). Defaults to None.
        """
        
        # Handle deprecated finite_elements parameter
        if finite_elements is not None:
            if finite_elements_number != 35 and finite_elements_number != finite_elements:
                # Check if finite_elements_number was explicitly set (not using default) and conflicts with finite_elements
                raise ValueError(FINITE_ELEMENTS_NUMBER_AND_FINITE_ELEMENTS_BOTH_SPECIFIED_ERROR)
            finite_elements_number = finite_elements
            print(FINITE_ELEMENTS_DEPRECATED_WARNING)
        
        # Use instance default if overwrite not explicitly provided
        if overwrite is None:
            overwrite = self.auto_confirm
        
        # Get default use_oep value if not provided
        if use_oep is None:
            use_oep = XC_FUNCTIONAL_OEP_DEFAULT.get(self.scf_xc_functional, False)
        
        # Safety confirmation prompt
        print("\n" + "="*75)
        print("WARNING: This script will generate/overwrite dataset files.")
        print("This operation may take a long time and will modify existing data in {}".format(self.data_root))
        print("="*75)
        if overwrite:
            print("\nAuto-confirming: Proceeding with data generation...")
        else:
            user_input = input("\nDo you want to proceed? (y/n): ").strip().lower()
            if user_input != 'y':
                print("Operation cancelled by user.")
                exit(0)
        print("\nStarting data generation...\n")

        self.generator.generate_data(
            # Required arguments
            data_root                   = self.data_root,
            atomic_number_list          = atomic_number_list,
            n_electrons_list            = n_electrons_list,
            use_oep                     = use_oep,
            scf_xc_functional           = self.scf_xc_functional,
            forward_pass_xc_functionals = self.forward_pass_xc_functional_list,

            # Arguments controlling the contents of the dataset
            save_energy_density         = save_energy_density,
            save_intermediate           = save_intermediate,
            save_full_spectrum          = save_full_spectrum,
            save_full_orbitals          = save_full_orbitals,
            save_derivative_matrix      = save_derivative_matrix,
            start_configuration_index   = start_configuration_index,

            # Arguments controlling the generation process
            # Grid, basis, and mesh parameters
            domain_size                 = domain_size,
            finite_elements_number      = finite_elements_number,
            polynomial_order            = polynomial_order,
            quadrature_point_number     = quadrature_point_number,
            oep_basis_number            = oep_basis_number,
            mesh_type                   = mesh_type,
            mesh_concentration          = mesh_concentration,
            mesh_spacing                = mesh_spacing,

            # SCF convergence parameters
            scf_tolerance               = scf_tolerance,
            max_scf_iterations          = max_scf_iterations,
            max_scf_iterations_outer    = max_scf_iterations_outer,
            use_pulay_mixing            = use_pulay_mixing,
            use_preconditioner          = use_preconditioner,
            pulay_mixing_parameter      = pulay_mixing_parameter,
            pulay_mixing_history        = pulay_mixing_history,
            pulay_mixing_frequency      = pulay_mixing_frequency,
            linear_mixing_alpha1        = linear_mixing_alpha1,
            linear_mixing_alpha2        = linear_mixing_alpha2,

            # Advanced functional parameters
            hybrid_mixing_parameter           = hybrid_mixing_parameter,
            frequency_quadrature_point_number = frequency_quadrature_point_number,
            angular_momentum_cutoff           = angular_momentum_cutoff,
            double_hybrid_flag                = double_hybrid_flag,
            oep_mixing_parameter              = oep_mixing_parameter,
            enable_parallelization            = enable_parallelization,

            # Debugging and verbose parameters
            verbose                     = verbose,
        )


    # Data loading methods
    def load_data(self,
        # Required arguments
        configuration_index_list    : Optional[List[int]] = None, # If None, load data for all configurations
        features_list               : List[str] = ["rho", "grad_rho", "lap_rho", "hartree", "lda_xc_potential"],

        # Control arguments
        use_radius_cutoff           : bool = False,
        use_feature_round_off       : bool = False,
        smooth_xc_data              : bool = False,
        close_shell_only            : bool = False,
        charge_neutral_only         : bool = False,
        spin_unpolarized_only       : bool = False,
        include_energy_density      : bool = False,
        include_intermediate        : bool = False,
        smooth_reduced_laplacian    : bool = False,
        smooth_rho_for_derivatives  : bool = False,
        print_debug_info            : bool = False,
        print_summary               : bool = False,

        # Additional arguments, for parameter 'tuning'
        radius_cutoff_rho_threshold : Optional[float] = None,
        radius_cutoff_v_x_threshold : Optional[float] = None,
        radius_cutoff_v_c_threshold : Optional[float] = None,
        smooth_radius_threshold     : Optional[float] = None,
        smooth_method               : Optional[str]   = None,
        smooth_kwargs               : Optional[Dict]  = None,

        # Log-linear asymptotic rho smoothing (see DataProcessor.smooth_rho_with_loglinear_asymptote; r_a < r_b < r_c in Bohr)
        smooth_rho_r_a              : Optional[float] = None,
        smooth_rho_r_b              : Optional[float] = None,
        smooth_rho_r_c              : Optional[float] = None,
        smooth_rho_log_floor        : Optional[float] = None,


        # Deprecated arguments
        use_cutoff                  : Literal[None] = None,        # Deprecated: use use_radius_cutoff instead
        cutoff_rho_threshold        : Literal[None] = None,        # Deprecated: use radius_cutoff_rho_threshold instead
        cutoff_v_x_threshold        : Literal[None] = None,        # Deprecated: use radius_cutoff_v_x_threshold instead
        cutoff_v_c_threshold        : Literal[None] = None,        # Deprecated: use radius_cutoff_v_c_threshold instead
        atomic_number_list          : Optional[List[int]] = None,  # Deprecated: use configuration_index_list instead
        smooth_vxc                  : Optional[bool] = None,       # Deprecated: use smooth_xc_data instead
    ) -> AtomicDataset:
        """
        Load training data for all configurations or specified ones.
        
        Parameters
        ----------
        configuration_index_list : Optional[List[int]]
            List of configuration indices (1-based) to load. If None, loads all available configurations.
            Atomic numbers are read from meta.json files in each atom folder.
        features_list : List[str]
            List of features to load
        
        use_radius_cutoff : bool
            Whether to apply cutoff filtering to truncate the radial grid
        use_feature_round_off : bool
            Whether to apply feature round-off (e.g., lower/upper bounds for small rho or potentials)
        smooth_xc_data : bool
            Whether to apply smoothing to V_xc and E_xc data
        close_shell_only : bool
            Whether to only load data for closed shell atoms
        charge_neutral_only : bool
            Whether to only load charge neutral configurations. For atom_XXX folders (old convention),
            defaults to charge neutral for backward compatibility.
        spin_unpolarized_only : bool
            Whether to only load configurations whose Aufbau-style occupations have
            ``occ_spin_up`` equal to ``occ_spin_down`` element-wise for the given atomic number
            and ``n_electrons`` read from each folder's metadata (same source as the data loader).
        include_energy_density : bool
            Whether to also load energy density data
        smooth_reduced_laplacian : bool
            Whether to smooth the reduced Laplacian at large radius using the weak-form 
        smooth_rho_for_derivatives : bool
            If True, apply log-linear asymptotic smoothing to rho and compute grad/lap (and s, q) from
            that rho via the FEM in out.txt instead of reading grad_rho.txt / lap_rho.txt.
        smooth_rho_r_a, smooth_rho_r_b, smooth_rho_r_c : Optional[float]
            Fit / blend radii in Bohr (``r_a < r_b < r_c`` for ``DataProcessor.smooth_rho_with_loglinear_asymptote``).
            If None, uses ``self.default_smooth_rho_r_a`` / ``self.default_smooth_rho_r_b`` / ``self.default_smooth_rho_r_c``
            (7.0, 10.0, 13.0 on a new manager).
        smooth_rho_log_floor : Optional[float]
            Floor for log fit in the asymptotic window. If None, uses ``self.default_smooth_rho_log_floor`` (1e-30 on a new manager).
        include_intermediate : bool
            Whether to also load data from intermediate iteration folders (outer_iter_XX)
        print_debug_info : bool
            Whether to print debug information
        
        radius_cutoff_rho_threshold : Optional[float]
            Threshold for rho data when use_radius_cutoff is True
        radius_cutoff_v_x_threshold : Optional[float]
            Threshold for v_x data when use_radius_cutoff is True
        radius_cutoff_v_c_threshold : Optional[float]
            Threshold for v_c data when use_radius_cutoff is True
        smooth_radius_threshold : Optional[float]
            Radius threshold for smoothing. Values with r > r_smooth_threshold will be smoothed.
            Default is 5.0.
        smooth_method : Optional[str]
            Smoothing method: 'lowpass' , 'savgol'(default), 'moving_avg', 'spline', 'gaussian', 'exp_weighted', 'cascade'
            - 'lowpass': Low-pass Butterworth filter (RECOMMENDED for high-frequency oscillations)
            - 'savgol': Savitzky-Golay filter (preserves data shape, but may not filter high-freq well)
            - 'moving_avg': Moving average (simple, good for high-frequency filtering)
            - 'spline': Spline smoothing (controllable smoothness)
            - 'gaussian': Gaussian filter (good smoothing, adjustable strength)
            - 'exp_weighted': Exponentially weighted moving average (smooth for large r)
            - 'cascade': Apply multiple smoothing methods in sequence (strongest filtering)
        smooth_kwargs : Optional[Dict]
            Additional parameters for smoothing methods:
            - lowpass: cutoff (default: 0.05), order (default: 6) - lower cutoff = stronger filtering
            - savgol: window_length (default: min(30% of data, len(data)//2*2+1)), polyorder (default: 2)
            - moving_avg: window_size (default: 25% of data, min 25) - larger = stronger filtering
            - spline: s (smoothing factor, default: len(data) * variance * 0.8) - larger = stronger
            - gaussian: sigma (default: max(2.0, 1% of data length)) - larger = stronger filtering
            - exp_weighted: alpha (default: 0.15) - smaller = stronger filtering
            - cascade: methods (list), kwargs_list (list of kwargs for each method)
    
        Deprecated Parameters
        --------------------
        atomic_number_list : Optional[List[int]]
            Deprecated: Use configuration_index_list instead. Atomic numbers are now read from meta.json files.
        smooth_vxc : Optional[bool]
            Deprecated: Use smooth_xc_data instead.

        Returns
        -------
        AtomicDataset: Atomic dataset
        """

        # handle deprecated arguments
        if use_cutoff is not None:
            warnings.warn(USE_CUTOFF_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            use_radius_cutoff = use_cutoff
        if cutoff_rho_threshold is not None:
            warnings.warn(CUTOFF_RHO_THRESHOLD_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            radius_cutoff_rho_threshold = cutoff_rho_threshold
        if cutoff_v_x_threshold is not None:
            warnings.warn(CUTOFF_V_X_THRESHOLD_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            radius_cutoff_v_x_threshold = cutoff_v_x_threshold
        if cutoff_v_c_threshold is not None:
            warnings.warn(CUTOFF_V_C_THRESHOLD_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            radius_cutoff_v_c_threshold = cutoff_v_c_threshold
        
        # Handle deprecated atomic_number_list parameter
        if atomic_number_list is not None:
            if configuration_index_list is not None:
                raise ValueError(ATOMIC_NUMBER_LIST_AND_CONFIGURATION_INDEX_LIST_BOTH_SPECIFIED_ERROR)
            print(ATOMIC_NUMBER_LIST_DEPRECATED_WARNING)
            configuration_index_list = atomic_number_list  # Convert to configuration_index_list

        # Handle deprecated smooth_vxc parameter
        if smooth_vxc is not None:
            warnings.warn(SMOOTH_VXC_DEPRECATED_WARNING, DeprecationWarning, stacklevel=2)
            smooth_xc_data = smooth_vxc

        # Find all configuration indices if not specified
        if configuration_index_list is None:
            configuration_index_list = []
            seen_indices = set()
            for item in os.listdir(self.data_root):
                item_path = os.path.join(self.data_root, item)
                if not os.path.isdir(item_path):
                    continue
                
                # Try new format first (configuration_XXX)
                if item.startswith('configuration_'):
                    try:
                        configuration_index = int(item.split('_')[1])
                        if configuration_index not in seen_indices:
                            configuration_index_list.append(configuration_index)
                            seen_indices.add(configuration_index)
                    except ValueError:
                        continue
                
                # Fallback to old format for backward compatibility (atom_XXX)
                elif item.startswith('atom_'):
                    try:
                        configuration_index = int(item.split('_')[1])
                        if configuration_index not in seen_indices:
                            configuration_index_list.append(configuration_index)
                            seen_indices.add(configuration_index)
                    except ValueError:
                        continue

            configuration_index_list = sorted(configuration_index_list)

        # set default parameters
        if radius_cutoff_rho_threshold is None:
            radius_cutoff_rho_threshold = self.default_radius_cutoff_rho_threshold
        if radius_cutoff_v_x_threshold is None:
            radius_cutoff_v_x_threshold = self.default_radius_cutoff_v_x_threshold
        if radius_cutoff_v_c_threshold is None:
            radius_cutoff_v_c_threshold = self.default_radius_cutoff_v_c_threshold
        if smooth_radius_threshold is None:
            smooth_radius_threshold = self.default_smooth_radius_threshold
        if smooth_method is None:
            smooth_method = self.default_smooth_method
        if smooth_kwargs is None:
            smooth_kwargs = self.default_smooth_kwargs

        if smooth_rho_r_a is None:
            smooth_rho_r_a = self.default_smooth_rho_r_a
        if smooth_rho_r_b is None:
            smooth_rho_r_b = self.default_smooth_rho_r_b
        if smooth_rho_r_c is None:
            smooth_rho_r_c = self.default_smooth_rho_r_c
        if smooth_rho_log_floor is None:
            smooth_rho_log_floor = self.default_smooth_rho_log_floor


        # Type and value check
        # Check configuration_index_list (list of integers >= 1)
        if not isinstance(configuration_index_list, list):
            raise TypeError("parameter 'configuration_index_list' must be a list, get {} instead.".format(type(configuration_index_list)))
        if not all(isinstance(idx, int) and idx >= 1 for idx in configuration_index_list):
            raise TypeError("parameter 'configuration_index_list' must be a list of integers >= 1, get {} instead.".format(configuration_index_list))
        
        AtomicDataManager.check_features_list(features_list)
        features_list = DataLoader.check_and_normalize_features_list(features_list)

        assert isinstance(use_radius_cutoff, bool), \
            USE_RADIUS_CUTOFF_NOT_BOOL_ERROR.format(type(use_radius_cutoff))
        assert isinstance(use_feature_round_off, bool), \
            USE_FEATURE_ROUND_OFF_NOT_BOOL_ERROR.format(type(use_feature_round_off))
        assert isinstance(smooth_xc_data, bool), \
            SMOOTH_XC_DATA_NOT_BOOL_ERROR.format(type(smooth_xc_data))
        assert isinstance(close_shell_only, bool), \
            CLOSE_SHELL_ONLY_NOT_BOOL_ERROR.format(type(close_shell_only))
        assert isinstance(spin_unpolarized_only, bool), \
            SPIN_UNPOLARIZED_ONLY_NOT_BOOL_ERROR.format(type(spin_unpolarized_only))
        assert isinstance(include_energy_density, bool), \
            INCLUDE_ENERGY_DENSITY_NOT_BOOL_ERROR.format(type(include_energy_density))
        assert isinstance(include_intermediate, bool), \
            INCLUDE_INTERMEDIATE_NOT_BOOL_ERROR.format(type(include_intermediate))
        assert isinstance(smooth_reduced_laplacian, bool), \
            SMOOTH_REDUCED_LAPLACIAN_NOT_BOOL_ERROR.format(type(smooth_reduced_laplacian))
        assert isinstance(smooth_rho_for_derivatives, bool), \
            SMOOTH_RHO_FOR_DERIVATIVES_NOT_BOOL_ERROR.format(type(smooth_rho_for_derivatives))
        assert isinstance(print_debug_info, bool), \
            PRINT_DEBUG_INFO_NOT_BOOL_ERROR.format(type(print_debug_info))
        assert isinstance(print_summary, bool), \
            PRINT_SUMMARY_NOT_BOOL_ERROR.format(type(print_summary))
    
        AtomicDataManager.check_radius_cutoff_threshold(
            radius_cutoff_rho_threshold,
            radius_cutoff_v_x_threshold,
            radius_cutoff_v_c_threshold,
        )
        AtomicDataManager.check_smooth_parameters(smooth_radius_threshold, smooth_method, smooth_kwargs)

        if smooth_reduced_laplacian and "lap_rho_reduced" not in features_list:
            warnings.warn(SMOOTH_REDUCED_LAPLACIAN_NO_LAP_RHO_REDUCED_WARNING, UserWarning, stacklevel=2)

        deriv_feats = {"grad_rho", "grad_rho_norm", "grad_rho_reduced", "lap_rho", "lap_rho_reduced"}
        if smooth_rho_for_derivatives and not deriv_feats.intersection(set(features_list)):
            warnings.warn(SMOOTH_RHO_FOR_DERIVATIVES_NO_DERIV_FEATURES_WARNING, UserWarning, stacklevel=2)
        
        # Filter to closed shell atoms only if requested
        if close_shell_only:
            filtered_config_indices = []
            for configuration_index in configuration_index_list:
                # Find configuration folder (with backward compatibility)
                config_folder = DataLoader._find_configuration_folder(self.data_root, configuration_index)
                if config_folder is None:
                    continue
                meta_path = os.path.join(config_folder, "meta.json")
                if os.path.exists(meta_path):
                    try:
                        import json
                        with open(meta_path, "r") as meta_file:
                            meta_data = json.load(meta_file)
                        z = int(meta_data.get("atomic_number", 0))
                        if z > 0:
                            occupation_info = OccupationInfo(z_nuclear=z, z_valence=z, all_electron_flag=True)
                            if occupation_info.closed_shell_flag:
                                filtered_config_indices.append(configuration_index)
                    except Exception as e:
                        print(f"Warning: Could not determine closed shell status for configuration {configuration_index}: {e}. Skipping...")
                        continue
            configuration_index_list = filtered_config_indices
            print(f"Filtered to {len(configuration_index_list)} closed shell configurations: {configuration_index_list}")
        

        # Filter to charge neutral atoms only if requested
        if charge_neutral_only:
            filtered_config_indices = []
            for configuration_index in configuration_index_list:
                # Find configuration folder (with backward compatibility)
                config_folder = DataLoader._find_configuration_folder(self.data_root, configuration_index)
                if config_folder is None:
                    continue
                folder_name = os.path.basename(config_folder)
                meta_path = os.path.join(config_folder, "meta.json")
                if os.path.exists(meta_path):
                    try:
                        import json
                        with open(meta_path, "r") as meta_file:
                            meta_data = json.load(meta_file)
                        z = int(meta_data.get("atomic_number", 0))
                        n_electrons = float(meta_data.get("n_electrons", 0))
                        if z > 0 and n_electrons == z:
                            filtered_config_indices.append(configuration_index)
                    except Exception as e:
                        print(f"Warning: Could not determine charge neutral status for configuration {configuration_index}: {e}. Skipping...")
                        continue

                # Backward compatibility: atom_XXX folders (old convention) default to charge neutral
                elif folder_name.startswith("atom_"):
                    filtered_config_indices.append(configuration_index)

            configuration_index_list = filtered_config_indices
            print(f"Filtered to {len(configuration_index_list)} charge neutral configurations: {configuration_index_list}")


        if spin_unpolarized_only:
            filtered_config_indices = []
            for configuration_index in configuration_index_list:
                config_folder = DataLoader._find_configuration_folder(self.data_root, configuration_index)
                if config_folder is None:
                    continue
                z, n_electrons, _ = DataLoader._load_metadata_from_meta(config_folder, configuration_index)
                if n_electrons is None:
                    print(
                        f"Warning: Could not determine n_electrons for spin-unpolarized filter "
                        f"for configuration {configuration_index}. Skipping..."
                    )
                    continue
                if z <= 0:
                    continue
                try:
                    occupation_info = OccupationInfo(
                        z_nuclear=z, z_valence=z, all_electron_flag=True, n_electrons=float(n_electrons)
                    )
                    if np.allclose(occupation_info.occ_spin_up, occupation_info.occ_spin_down):
                        filtered_config_indices.append(configuration_index)
                except Exception as e:
                    print(
                        f"Warning: Could not determine spin-unpolarized status for configuration "
                        f"{configuration_index}: {e}. Skipping..."
                    )
                    continue
            configuration_index_list = filtered_config_indices
            print(
                f"Filtered to {len(configuration_index_list)} spin-unpolarized configurations: "
                f"{configuration_index_list}"
            )


        # Load data for each configuration
        configuration_data_list, skipped_atoms = DataLoader.load_data(
            data_root                   = self.data_root,
            scf_xc_functional           = self.scf_xc_functional,
            forward_pass_xc_functionals = self.forward_pass_xc_functional_list,
            features_list               = features_list,
            configuration_index_list    = configuration_index_list,
            use_radius_cutoff           = use_radius_cutoff,
            use_feature_round_off       = use_feature_round_off,
            include_energy_density      = include_energy_density,
            include_intermediate        = include_intermediate,
            print_debug_info            = print_debug_info,
            radius_cutoff_rho_threshold = radius_cutoff_rho_threshold,
            radius_cutoff_v_x_threshold = radius_cutoff_v_x_threshold,
            radius_cutoff_v_c_threshold = radius_cutoff_v_c_threshold,
            smooth_reduced_laplacian    = smooth_reduced_laplacian,
            smooth_rho_for_derivatives  = smooth_rho_for_derivatives,
            smooth_rho_r_a              = smooth_rho_r_a,
            smooth_rho_r_b              = smooth_rho_r_b,
            smooth_rho_r_c              = smooth_rho_r_c,
            smooth_rho_log_floor        = smooth_rho_log_floor,
        )

        # Load shared derivative matrix if any configuration uses it
        shared_derivative_matrix = None
        shared_derivative_matrix_path = os.path.join(self.data_root, "derivative_matrix.npy")
        if os.path.exists(shared_derivative_matrix_path):
            # Check if any configuration uses shared derivative matrix
            for config_data in configuration_data_list:
                if hasattr(config_data, 'derivative_matrix_use_shared') and config_data.derivative_matrix_use_shared:
                    # Load shared derivative matrix once
                    shared_derivative_matrix = np.load(shared_derivative_matrix_path)
                    break

        dataset = AtomicDataset(
            # Basic attributes
            data_root                       = self.data_root,
            scf_xc_functional               = self.scf_xc_functional,
            forward_pass_xc_functional_list = self.forward_pass_xc_functional_list,
            features_list                   = features_list,
            
            # Other attributes
            radius_cutoff_rho_threshold     = radius_cutoff_rho_threshold,
            radius_cutoff_v_x_threshold     = radius_cutoff_v_x_threshold,
            radius_cutoff_v_c_threshold     = radius_cutoff_v_c_threshold,

            # Data attributes
            configuration_data_list         = configuration_data_list,
            shared_derivative_matrix        = shared_derivative_matrix,
        )

        if smooth_xc_data:
            dataset.smooth_xc_data(smooth_radius_threshold, smooth_method, smooth_kwargs)

        if print_debug_info or print_summary:
            dataset.print_info()
        
        return dataset



    def inverse_transform_features(
        self,
        X_transformed    : np.ndarray,
        scaler_X         : Optional[Any]  = None,
        transform_params : Optional[Dict] = None,
        feature_idx      : Optional[int]  = None
    ) -> np.ndarray:
        """
        Inverse transform features back to physical space.
        """
        return DataProcessor.inverse_transform_features(X_transformed, scaler_X, transform_params, feature_idx)


    def inverse_transform_predictions(
        self,
        y_pred           : np.ndarray,
        scaler_y         : Optional[Any] = None,
        transform_params : Optional[Dict] = None
    ) -> np.ndarray:
        """
        Inverse transform predictions back to physical space.
        """
        return DataProcessor.inverse_transform_predictions(y_pred, scaler_y, transform_params)


    # Type and value checks
    @staticmethod
    def check_data_root(data_root) -> None:
        """
        Check if the data root is a valid directory.
        """
        if not isinstance(data_root, str):
            raise TypeError(DATA_ROOT_NOT_STRING_ERROR.format(type(data_root)))
        if not os.path.exists(data_root):
            raise FileNotFoundError(DATA_ROOT_NOT_EXIST_ERROR.format(data_root))


    @staticmethod
    def check_scf_xc_functional(scf_xc_functional) -> None:
        """
        Check if the SCF XC functional is a valid string.
        """
        from ..xc.functional_requirements import VALID_XC_FUNCTIONAL_LIST
        if not isinstance(scf_xc_functional, str):
            raise TypeError(SCF_XC_FUNCTIONAL_NOT_STRING_ERROR.format(type(scf_xc_functional)))
        if scf_xc_functional not in VALID_XC_FUNCTIONAL_LIST:
            raise ValueError(SCF_XC_FUNCTIONAL_NOT_IN_VALID_LIST_ERROR.format(VALID_XC_FUNCTIONAL_LIST, scf_xc_functional))


    @staticmethod
    def check_forward_pass_xc_functional_list(forward_pass_xc_functional_list) -> None:
        """
        Check if the forward pass XC functional list is a valid list of strings.
        """
        from ..xc.functional_requirements import VALID_XC_FUNCTIONAL_LIST
        if not isinstance(forward_pass_xc_functional_list, list):
            raise TypeError(FORWARD_PASS_XC_FUNCTIONAL_LIST_NOT_LIST_ERROR.format(type(forward_pass_xc_functional_list)))
        if not all(isinstance(functional, str) for functional in forward_pass_xc_functional_list):
            raise TypeError(FORWARD_PASS_XC_FUNCTIONAL_LIST_NOT_LIST_OF_STRINGS_ERROR.format(type(forward_pass_xc_functional_list)))
        if not all(functional in VALID_XC_FUNCTIONAL_LIST for functional in forward_pass_xc_functional_list):
            raise ValueError(FORWARD_PASS_XC_FUNCTIONAL_LIST_NOT_IN_VALID_LIST_ERROR.format(VALID_XC_FUNCTIONAL_LIST, forward_pass_xc_functional_list))


    @staticmethod
    def check_features_list(features_list) -> None:
        """
        Check if the features list is a valid list of strings.
        """
        if not isinstance(features_list, list):
            raise TypeError(FEATURES_LIST_NOT_LIST_ERROR.format(type(features_list)))
        if not all(isinstance(feature, str) for feature in features_list):
            raise TypeError(FEATURES_LIST_NOT_LIST_OF_STRINGS_ERROR.format(type(features_list)))
        for feature in features_list:
            normalized_feature = FEATURE_ALIASES.get(feature, feature)
            if normalized_feature not in NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL:
                raise ValueError(format_invalid_feature_error(feature, NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL))



    @staticmethod
    def check_radius_cutoff_threshold(
        radius_cutoff_rho_threshold,
        radius_cutoff_v_x_threshold,
        radius_cutoff_v_c_threshold,
    ) -> None:
        """
        Check if the cutoff rho threshold is a valid float.
        """

        # Cutoff rho threshold
        if not isinstance(radius_cutoff_rho_threshold, float):
            raise TypeError(RADIUS_CUTOFF_RHO_THRESHOLD_NOT_FLOAT_ERROR.format(type(radius_cutoff_rho_threshold)))
        if radius_cutoff_rho_threshold <= 0:
            raise ValueError(RADIUS_CUTOFF_RHO_THRESHOLD_NOT_POSITIVE_ERROR.format(radius_cutoff_rho_threshold))
        
        # Cutoff v_x threshold
        if radius_cutoff_v_x_threshold is not None:
            if not isinstance(radius_cutoff_v_x_threshold, float):
                raise TypeError(RADIUS_CUTOFF_V_X_THRESHOLD_NOT_FLOAT_ERROR.format(type(radius_cutoff_v_x_threshold)))
            if radius_cutoff_v_x_threshold <= 0:
                raise ValueError(RADIUS_CUTOFF_V_X_THRESHOLD_NOT_POSITIVE_ERROR.format(radius_cutoff_v_x_threshold))
        
        # Cutoff v_c threshold
        if radius_cutoff_v_c_threshold is not None:
            if not isinstance(radius_cutoff_v_c_threshold, float):
                raise TypeError(RADIUS_CUTOFF_V_C_THRESHOLD_NOT_FLOAT_ERROR.format(type(radius_cutoff_v_c_threshold)))
            if radius_cutoff_v_c_threshold <= 0:
                raise ValueError(RADIUS_CUTOFF_V_C_THRESHOLD_NOT_POSITIVE_ERROR.format(radius_cutoff_v_c_threshold))

    @staticmethod
    def check_smooth_parameters(smooth_radius_threshold, smooth_method, smooth_kwargs) -> None:
        """
        Check if the smooth radius threshold is a valid float.
        """
        # Smooth radius threshold
        if not isinstance(smooth_radius_threshold, float):
            raise TypeError(SMOOTH_RADIUS_THRESHOLD_NOT_FLOAT_ERROR.format(type(smooth_radius_threshold)))
        if smooth_radius_threshold <= 0:
            raise ValueError(SMOOTH_RADIUS_THRESHOLD_NOT_POSITIVE_ERROR.format(smooth_radius_threshold))

        # Smooth method
        if not isinstance(smooth_method, str):
            raise TypeError(SMOOTH_METHOD_NOT_STRING_ERROR.format(type(smooth_method)))
        if smooth_method not in VALID_SMOOTH_METHODS:
            raise ValueError(SMOOTH_METHOD_NOT_IN_VALID_LIST_ERROR.format(VALID_SMOOTH_METHODS, smooth_method))

        # Smooth kwargs
        if not isinstance(smooth_kwargs, dict):
            raise TypeError(SMOOTH_KWARGS_NOT_DICT_ERROR.format(type(smooth_kwargs)))


    @staticmethod
    def check_configuration_data_list(configuration_data_list) -> None:
        """
        Check if the configuration data list is a valid list of SingleConfigurationData.
        Uses class name check to avoid importlib.reload issues (isinstance fails across reloads).
        """
        if not isinstance(configuration_data_list, list):
            raise TypeError(CONFIGURATION_DATA_LIST_NOT_LIST_ERROR.format(type(configuration_data_list)))
        def _is_single_config(cfg) -> bool:
            return type(cfg).__name__ == "SingleConfigurationData"
        if not all(_is_single_config(c) for c in configuration_data_list):
            invalid_types = {
                type(c) for c in configuration_data_list
                if not _is_single_config(c)
            }
            raise TypeError(
                CONFIGURATION_DATA_LIST_NOT_LIST_OF_SINGLE_CONFIGURATION_DATA_ERROR.format(invalid_types)
            )


    @staticmethod
    def check_atomic_number_list(atomic_number_list: List[int]) -> None:
        """
        Check if the atomic number list is a valid list of integers.
        """
        if not isinstance(atomic_number_list, list):
            raise TypeError(ATOMIC_NUMBER_LIST_NOT_LIST_ERROR.format(type(atomic_number_list)))
        if not all(isinstance(atomic_number, int) for atomic_number in atomic_number_list):
            raise TypeError(ATOMIC_NUMBER_LIST_NOT_LIST_OF_INTEGERS_ERROR.format(type(atomic_number_list)))
        if not all(atomic_number >= 1 and atomic_number <= 92 for atomic_number in atomic_number_list):
            raise ValueError(ATOMIC_NUMBER_LIST_NOT_IN_VALID_RANGE_ERROR.format(atomic_number_list))


    # Default parameter 'setting' methods
    def set_default_radius_cutoff_rho_threshold(self, radius_cutoff_rho_threshold: float):
        assert isinstance(radius_cutoff_rho_threshold, float), \
            RADIUS_CUTOFF_RHO_THRESHOLD_NOT_FLOAT_ERROR.format(type(radius_cutoff_rho_threshold))
        assert radius_cutoff_rho_threshold > 0, \
            RADIUS_CUTOFF_RHO_THRESHOLD_NOT_POSITIVE_ERROR.format(radius_cutoff_rho_threshold)
        self.default_radius_cutoff_rho_threshold = radius_cutoff_rho_threshold

    def set_default_radius_cutoff_v_x_threshold(self, radius_cutoff_v_x_threshold: float):
        assert isinstance(radius_cutoff_v_x_threshold, float), \
            RADIUS_CUTOFF_V_X_THRESHOLD_NOT_FLOAT_ERROR.format(type(radius_cutoff_v_x_threshold))
        assert radius_cutoff_v_x_threshold > 0, \
            RADIUS_CUTOFF_V_X_THRESHOLD_NOT_POSITIVE_ERROR.format(radius_cutoff_v_x_threshold)
        self.default_radius_cutoff_v_x_threshold = radius_cutoff_v_x_threshold

    def set_default_radius_cutoff_v_c_threshold(self, radius_cutoff_v_c_threshold: float):
        assert isinstance(radius_cutoff_v_c_threshold, float), \
            RADIUS_CUTOFF_V_C_THRESHOLD_NOT_FLOAT_ERROR.format(type(radius_cutoff_v_c_threshold))
        assert radius_cutoff_v_c_threshold > 0, \
            RADIUS_CUTOFF_V_C_THRESHOLD_NOT_POSITIVE_ERROR.format(radius_cutoff_v_c_threshold)
        self.default_radius_cutoff_v_c_threshold = radius_cutoff_v_c_threshold

    def set_default_smooth_radius_threshold(self, smooth_radius_threshold: float):
        assert isinstance(smooth_radius_threshold, float), \
            SMOOTH_RADIUS_THRESHOLD_NOT_FLOAT_ERROR.format(type(smooth_radius_threshold))
        assert smooth_radius_threshold > 0, \
            SMOOTH_RADIUS_THRESHOLD_NOT_POSITIVE_ERROR.format(smooth_radius_threshold)
        self.default_smooth_radius_threshold = smooth_radius_threshold
    
    def set_default_smooth_method(self, smooth_method: str):
        assert isinstance(smooth_method, str), \
            SMOOTH_METHOD_NOT_STRING_ERROR.format(type(smooth_method))
        assert smooth_method in VALID_SMOOTH_METHODS, \
            SMOOTH_METHOD_NOT_IN_VALID_LIST_ERROR.format(VALID_SMOOTH_METHODS, smooth_method)
        self.default_smooth_method = smooth_method
    
    def set_default_smooth_kwargs(self, smooth_kwargs: Dict):
        assert isinstance(smooth_kwargs, dict), \
            SMOOTH_KWARGS_NOT_DICT_ERROR.format(type(smooth_kwargs))
        self.default_smooth_kwargs = smooth_kwargs
    

    @staticmethod
    def integrate_energy_density(r: np.ndarray, w: np.ndarray, energy_density: np.ndarray) -> float:
        """
        Integrate energy density over the grid.
        """
        return AtomicDataset.integrate_energy_density(r, w, energy_density)