
__author__ = "Qihao Cheng"


"""Data loading utilities for atomic DFT calculations."""
import numpy as np
import os
import glob
import json
import warnings
from dataclasses import dataclass
from typing import Tuple, List, Optional, Dict, Literal, Union, get_args, TYPE_CHECKING

from ..xc.lda import lda_exchange_potential_generic, lda_correlation_potential_generic
from .data_processing import DataProcessor

if TYPE_CHECKING:
    from ..solver import AtomicDFTSolver
# format_error_message is imported locally inside functions to avoid circular import



# Error messages for loading data
REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR = \
    "parameter 'reference_array_length' must be an integer, get {} instead."
SCF_FOLDER_NOT_FOUND_ERROR = \
    "SCF folder {} does not exist."
SCF_FOLDER_PATH_NOT_STRING_ERROR = \
    "SCF folder path must be a string, get {} instead."

FORWARD_PASS_FOLDER_NOT_FOUND_ERROR = \
    "Forward pass folder {} does not exist."
FORWARD_PASS_FOLDER_PATH_NOT_STRING_ERROR = \
    "Forward pass folder path must be a string, get {} instead."
FORWARD_PASS_FOLDER_PATH_LIST_NOT_LIST_ERROR = \
    "Forward pass folder path list must be a list, get {} instead."

QUADRATURE_NODES_NOT_NUMPY_ARRAY_ERROR = \
    "Quadrature nodes must be a numpy array, get {} instead."
QUADRATURE_NODES_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "Quadrature nodes length {} does not match reference array length {}."
QUADRATURE_NODES_NOT_1D_ARRAY_ERROR = \
    "Quadrature nodes must be a 1D numpy array, get {} instead."

QUADRATURE_WEIGHTS_NOT_NUMPY_ARRAY_ERROR = \
    "Quadrature weights must be a numpy array, get {} instead."
QUADRATURE_WEIGHTS_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "Quadrature weights length {} does not match reference array length {}."
QUADRATURE_WEIGHTS_NOT_1D_ARRAY_ERROR = \
    "Quadrature weights must be a 1D numpy array, get {} instead."

RHO_DATA_MUST_BE_INCLUDED_IN_FEATURES_LIST_ERROR = \
    "Rho data must be included in the provided features list."
RHO_NOT_NUMPY_ARRAY_ERROR = \
    "Rho must be a numpy array, get {} instead."
RHO_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "Rho length {} does not match reference array length {}."
RHO_NOT_1D_ARRAY_ERROR = \
    "Rho must be a 1D numpy array, get {} instead."

GRAD_RHO_NOT_NUMPY_ARRAY_ERROR = \
    "The gradient of rho must be a numpy array, get {} instead."
GRAD_RHO_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "The gradient of rho's length {} does not match reference array length {}."
GRAD_RHO_NOT_1D_ARRAY_ERROR = \
    "The gradient of rho must be a 1D numpy array, get {} instead."

LAP_RHO_NOT_NUMPY_ARRAY_ERROR = \
    "The Laplacian of rho must be a numpy array, get {} instead."
LAP_RHO_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "The Laplacian of rho's length {} does not match reference array length {}."
LAP_RHO_NOT_1D_ARRAY_ERROR = \
    "The Laplacian of rho must be a 1D numpy array, get {} instead."

HARTREE_NOT_NUMPY_ARRAY_ERROR = \
    "Hartree must be a numpy array, get {} instead."
HARTREE_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "Hartree length {} does not match reference array length {}."
HARTREE_NOT_1D_ARRAY_ERROR = \
    "Hartree must be a 1D numpy array, get {} instead."

V_X_NOT_NUMPY_ARRAY_ERROR = \
    "V_x must be a numpy array, get {} instead."
V_X_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "V_x length {} does not match reference array length {}."
V_X_NOT_1D_ARRAY_ERROR = \
    "V_x must be a 1D numpy array, get {} instead."

V_C_NOT_NUMPY_ARRAY_ERROR = \
    "V_c must be a numpy array, get {} instead."
V_C_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "V_c length {} does not match reference array length {}."
V_C_NOT_1D_ARRAY_ERROR = \
    "V_c must be a 1D numpy array, get {} instead."

E_X_NOT_NUMPY_ARRAY_ERROR = \
    "E_x must be a numpy array, get {} instead."
E_X_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "E_x length {} does not match reference array length {}."
E_X_NOT_1D_ARRAY_ERROR = \
    "E_x must be a 1D numpy array, get {} instead."

E_C_NOT_NUMPY_ARRAY_ERROR = \
    "E_c must be a numpy array, get {} instead."
E_C_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "E_c length {} does not match reference array length {}."
E_C_NOT_1D_ARRAY_ERROR = \
    "E_c must be a 1D numpy array, get {} instead."

FEATURES_DATA_LENGTH_NOT_EQUAL_TO_FEATURES_LIST_LENGTH_ERROR = \
    "Features data length {} does not match features list length {}."
FEATURES_DATA_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR = \
    "Features data length {} does not match reference array length {}."

CUTOFF_INDEX_NOT_FOUND_ERROR = \
    "Cutoff index not found, please check the filtering criteria."
V_X_FILE_NOT_FOUND_ERROR = \
    "Could not find v_x.txt or any matching v_x_*.txt file (excluding *_uniform.txt) in {}."
V_C_FILE_NOT_FOUND_ERROR = \
    "Could not find v_c.txt or any matching v_c_*.txt file (excluding *_uniform.txt) in {}."

E_X_FILE_NOT_FOUND_ERROR = \
    "Could not find e_x.txt or any matching e_x_*.txt file (excluding *_uniform.txt) in {}."
E_C_FILE_NOT_FOUND_ERROR = \
    "Could not find e_c.txt or any matching e_c_*.txt file (excluding *_uniform.txt) in {}."


# Error messages for SingleConfigurationData
FEATURES_DATA_NOT_NUMPY_ARRAY_ERROR = \
    "Features must be a numpy array, get {} instead."
SCF_XC_DATA_NOT_TUPLE_ERROR = \
    "SCF XC data must be a tuple, get {} instead."
FORWARD_PASS_XC_DATA_LIST_NOT_LIST_ERROR = \
    "Forward pass XC data list must be a list, get {} instead."
ATOMIC_NUMBERS_NOT_NUMPY_ARRAY_ERROR = \
    "Atomic numbers must be a numpy array, get {} instead."
QUADRATURE_POINTS_NOT_NUMPY_ARRAY_ERROR = \
    "Quadrature points must be a numpy array, get {} instead."
CONFIGURATION_IDS_NOT_NUMPY_ARRAY_ERROR = \
    "Configuration IDs must be a numpy array, get {} instead."
QUADRATURE_NODES_FILTERED_NOT_NUMPY_ARRAY_ERROR = \
    "Filtered quadrature nodes must be a numpy array, get {} instead."
CUTOFF_RADIUS_NOT_FLOAT_ERROR = \
    "Cutoff radius must be a float, get {} instead."
CUTOFF_IDX_NOT_INT_ERROR = \
    "Cutoff index must be an integer, get {} instead."
CONFIGURATION_ID_NOT_INT_ERROR = \
    "Configuration ID must be an integer, get {} instead."
N_POINTS_NOT_INT_ERROR = \
    "Number of points must be an integer, get {} instead."
N_FILTERED_NOT_INT_ERROR = \
    "Number of filtered points must be an integer, get {} instead."
FOLDER_LABEL_NOT_STRING_ERROR = \
    "Folder label must be a string, get {} instead."

ATOMIC_NUMBER_NOT_INTEGER_ERROR = \
    "Atomic number must be an integer, get {} instead."
N_ELECTRONS_NOT_FLOAT_ERROR = \
    "Number of electrons must be a float, get {} instead."
FEATURES_LIST_NOT_LIST_ERROR = \
    "Features list must be a list, get {} instead."
FEATURES_DATA_NOT_2D_ARRAY_ERROR = \
    "Features data must be a 2D numpy array, get {} instead."
SCF_XC_DATA_NOT_1D_ARRAY_ERROR = \
    "SCF XC data must be a 1D numpy array, get {} instead."
FORWARD_PASS_XC_DATA_NOT_1D_ARRAY_ERROR = \
    "Forward pass XC data must be a 1D numpy array, get {} instead."
ATOMIC_NUMBERS_NOT_1D_ARRAY_ERROR = \
    "Atomic numbers must be a 1D numpy array, get {} instead."
QUADRATURE_POINTS_NOT_1D_ARRAY_ERROR = \
    "Quadrature points must be a 1D numpy array, get {} instead."
CONFIGURATION_IDS_NOT_1D_ARRAY_ERROR = \
    "Configuration IDs must be a 1D numpy array, get {} instead."
QUADRATURE_NODES_FILTERED_NOT_1D_ARRAY_ERROR = \
    "Filtered quadrature nodes must be a 1D numpy array, get {} instead."
CUTOFF_RADIUS_NOT_POSITIVE_ERROR = \
    "Cutoff radius must be positive, get {} instead."
CUTOFF_IDX_NOT_NON_NEGATIVE_ERROR = \
    "Cutoff index must be non-negative, get {} instead."
N_POINTS_NOT_POSITIVE_ERROR = \
    "Number of points must be positive, get {} instead."
N_FILTERED_NOT_POSITIVE_ERROR = \
    "Number of filtered points must be positive, get {} instead."

WEIGHTS_DATA_NOT_NUMPY_ARRAY_ERROR = \
    "parameter 'weights_data' must be a numpy array, get {} instead."
WEIGHTS_DATA_NOT_1D_OR_2D_ARRAY_ERROR = \
    "parameter 'weights_data' must be a 1D or 2D numpy array, get {} instead."
WEIGHTS_DATA_LENGTH_NOT_EQUAL_TO_N_FILTERED_ERROR = \
    "parameter 'weights_data' length {} does not match number of filtered points {}."

NO_DATA_LOADED_ERROR = \
    "No data was successfully loaded!."

FUNCTIONAL_NOT_IN_CONFIGURATION_DATA_ERROR = \
    "Functional '{}' is not in configuration data: scf={}, forward_pass={}."
TARGET_COMPONENT_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'target_component' must be in {}, get {} instead."


# Warning messages
SCF_FOLDER_NOT_FOUND_WARNING = \
    "WARNING: SCF folder {} does not exist, skipping atom {}..."
FORWARD_PASS_FOLDER_NOT_FOUND_WARNING = \
    "WARNING: Forward pass folder {} does not exist, skipping atom {}..."
CONVERGED_CONFIGURATION_DATA_NOT_FOUND_WARNING = \
    "WARNING: Converged configuration data not found for atom {}, skipping..."
INTERMEDIATE_CONFIGURATION_DATA_NOT_FOUND_WARNING = \
    "WARNING: Intermediate configuration data not found for atom {} in folder {}, skipping..."
CONFIGURATION_FOLDER_NOT_FOUND_WARNING = \
    "WARNING: Configuration folder not found for index {0} (tried both configuration_{0:03d} and atom_{0:03d})"
ATOMIC_NUMBER_EXTRACTION_FAILED_WARNING = \
    "WARNING: Cannot extract atomic_number from folder name '{}' (expected format: atom_XXX). Using configuration index {} as fallback."
META_JSON_NOT_FOUND_WARNING = \
    "WARNING: meta.json not found in folder '{}' and folder name doesn't match 'atom_XXX' format. Using configuration index {} as fallback."
CONFIGURATION_NOT_CONVERGED_WARNING = \
    "WARNING: Configuration {} (atomic_number={}) did not converge, skipping..."
SKIPPED_ATOMS_WARNING = \
    "WARNING: Skipped {} atom(s) due to missing data: {}"
ATOMIC_NUMBER_NOT_FOUND_IN_META_WARNING = \
    "WARNING: meta.json not found in {} or atomic_number is invalid. Using configuration index {} as temporary atomic number."
SMOOTH_RHO_NON_MONOTONE_DROPPED_WARNING = \
    "WARNING: smooth_rho_for_derivatives=True but rho is not monotone non-increasing on the radial grid " \
    "(configuration_index={}, folder={}, atomic_number={}). Dropping this configuration."

POTENTIAL_WEIGHTS_DATA_ALREADY_SET_WARNING = \
    "WARNING: Potential weights data already set for atom {}, configuration id {}, update again will overwrite the existing weights data."
ENERGY_WEIGHTS_DATA_ALREADY_SET_WARNING = \
    "WARNING: Energy weights data already set for atom {}, configuration id {}, update again will overwrite the existing weights data."

# Debug messages
CONVERGED_CONFIGURATION_DATA_LOADED_DEBUG = \
    "Loaded converged configuration data [index: {:>3d}, atomic_number: {:>3d}, n_electrons: {:.6f}]: {:>5d} -> {:>5d} grid points (cutoff radius: {:.6f})"
INTERMEDIATE_CONFIGURATION_DATA_LOADED_DEBUG = \
    "\t Loaded intermediate configuration data for atom {:>3d}: {:>5d} -> {:>5d} grid points (number of electrons: {:.6f}, cutoff radius: {:.6f})"
FOUND_INTERMEDIATE_CONFIGURATION_DATA_DEBUG = \
    "Found {:>2d} intermediate configuration data for atom {:>3d}"



# Valid features list
NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL: List[str] = [
    "rho",              # electron density
    "grad_rho",         # density gradient
    "grad_rho_norm",    # density gradient's magnitude
    "grad_rho_reduced", # reduced gradient
    "lap_rho",          # Laplacian of density
    "lap_rho_reduced",  # reduced Laplacian
    "hartree",          # Hartree potential
    "lda_xc_potential", # LDA XC potential
]

NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY: List[str] = [
    "rho",
    "grad_rho",
    "grad_rho_norm",
    "grad_rho_reduced",
    "lap_rho",
    "lap_rho_reduced",
]

# Feature aliases that map to the valid feature list.
FEATURE_ALIASES: Dict[str, str] = {
    "s": "grad_rho_reduced",
    "q": "lap_rho_reduced",
    "grad_rho_abs": "grad_rho_norm",
    "grad_rho_mag": "grad_rho_norm",
    "lda_xc": "lda_xc_potential",  # backward compatibility
}

FEATURE_DESCRIPTIONS_DICT: Dict[str, str] = {
    "rho"              : "Electron density.",
    "grad_rho"         : "Density gradient.",
    "grad_rho_norm"    : "Magnitude of the density gradient.",
    "grad_rho_reduced" : "Reduced density gradient.",
    "lap_rho"          : "Laplacian of the density.",
    "lap_rho_reduced"  : "Reduced Laplacian of the density.",
    "hartree"          : "Hartree potential.",
    "lda_xc_potential" : "LDA exchange-correlation potential.",
    "lda_xc"           : "Alias of 'lda_xc_potential'.",
    "s"                : "Alias of 'grad_rho_reduced'.",
    "q"                : "Alias of 'lap_rho_reduced'.",
    "grad_rho_abs"     : "Alias of 'grad_rho_norm'.",
    "grad_rho_mag"     : "Alias of 'grad_rho_norm'.",
}



def format_invalid_feature_error(
    feature: str,
    valid_features_list: Optional[List[str]] = None,
) -> str:
    """
    Format an invalid feature error with a standardized description list.
    """
    if valid_features_list is None:
        keys = list(FEATURE_DESCRIPTIONS_DICT.keys())
    else:
        keys = list(valid_features_list)
        for alias, target in FEATURE_ALIASES.items():
            if target in valid_features_list and alias not in keys:
                keys.append(alias)
    max_key_len = max(len(key) for key in keys)
    descriptions = "\n".join(
        f"\t - {key:<{max_key_len}} : {FEATURE_DESCRIPTIONS_DICT[key]}"
        for key in keys
        if key in FEATURE_DESCRIPTIONS_DICT
    )
    return (
        f"Invalid feature: '{feature}', must be one of the following: \n"
        f"{descriptions}"
    )


def _expand_feature_aliases(valid_features: List[str]) -> List[str]:
    expanded = list(valid_features)
    for alias, target in FEATURE_ALIASES.items():
        if target in valid_features and alias not in expanded:
            expanded.append(alias)
    return expanded


# Full valid list (includes aliases)
VALID_FEATURES_LIST_FOR_POTENTIAL: List[str] = _expand_feature_aliases(
    NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL
)
VALID_FEATURES_LIST_FOR_ENERGY_DENSITY: List[str] = _expand_feature_aliases(
    NORMALIZED_VALID_FEATURES_LIST_FOR_ENERGY_DENSITY
)



# Lower / upper bounds for features data
RHO_LOWER_BOUND      : float =  1e-13   # In principle, rho should be positive, so we only consider the lower bound
GRAD_RHO_UPPER_BOUND : float = -1e-13   # In principle, grad_rho should be negative, so we only consider the upper bound
HARTREE_LOWER_BOUND  : float =  1e-13   # In principle, hartree potential should be positive, so we only consider the lower bound
LDA_XC_POTENTIAL_UPPER_BOUND : float = -1e-13   # In principle, lda_xc_potential should be positive, so we only consider the lower bound


# Valid target components for potential and energy
PotentialTargetComponent = Literal["v_xc", "v_x", "v_c", "v_x_v_c", "enhancement_factor_x"]
EnergyTargetComponent    = Literal["e_xc", "e_x", "e_c", "e_x_e_c", "enhancement_factor_x"]

VALID_POTENTIAL_TARGET_COMPONENTS = set(get_args(PotentialTargetComponent))
VALID_ENERGY_TARGET_COMPONENTS    = set(get_args(EnergyTargetComponent))


# Type aliases
# (v_x, v_c, e_x, e_c), where e_x and e_c are only available if include_energy_density is True
XCDataType = Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]



@dataclass
class SingleConfigurationData:
    """
    Data class for single configuration data.
    
    Attributes:
        atomic_number                   : int | float
        n_electrons                     : int | float
        configuration_id                : int
        cutoff_radius                   : float
        cutoff_idx                      : int
        n_points                        : int
        n_filtered                      : int
        folder_label                    : str
        scf_xc_functional               : str
        scf_folder_path                 : str
        forward_pass_xc_functional_list : List[str]
        forward_pass_folder_path_list   : List[str]
        features_list                   : List[str]

        features_data                   : np.ndarray of shape (n_filtered, n_features)
        quadrature_nodes                : np.ndarray of shape (n_points, )
        quadrature_weights              : np.ndarray of shape (n_points, )
        scf_xc_data                     : XCDataType (v_x, v_c, e_x, e_c), each of shape (n_filtered, )
        forward_pass_xc_data_list       : List[XCDataType], each of shape (n_filtered, )
        derivative_matrix_use_shared    : bool
            Whether this atom uses shared derivative matrix. If True, derivative_matrix is None and shared matrix is in dataset.
        derivative_matrix               : Optional[np.ndarray] of shape (n_elements, n_quad, n_quad)
            Derivative matrix (None if use_shared=True, otherwise stores unique matrix for this atom).
        potential_weights_data          : Optional[np.ndarray] of shape (n_filtered, n_targets)
        energy_weights_data             : Optional[np.ndarray] of shape (n_filtered, n_targets)
    """
    atomic_number                   : int | float
    n_electrons                     : int | float
    configuration_id                : int
    cutoff_radius                   : float
    cutoff_idx                      : int
    n_points                        : int
    n_filtered                      : int
    folder_label                    : str
    scf_xc_functional               : str
    scf_folder_path                 : str
    forward_pass_xc_functional_list : List[str]
    forward_pass_folder_path_list   : List[str]
    features_list                   : List[str]

    features_data                   : np.ndarray
    quadrature_nodes                : np.ndarray
    quadrature_weights              : np.ndarray
    scf_xc_data                     : XCDataType
    forward_pass_xc_data_list       : List[XCDataType]
    derivative_matrix_use_shared    : bool = False  # Whether this atom uses shared derivative matrix
    derivative_matrix               : Optional[np.ndarray] = None  # Derivative matrix (None if use_shared=True, otherwise stores unique matrix)


    def __post_init__(self):
        self._check_inputs_type_value_and_shape()
        
        # Initialize weights data
        self.potential_weights_data : Optional[np.ndarray] = None
        self.energy_weights_data    : Optional[np.ndarray] = None



    def set_potential_weights_data(self, weights_data: np.ndarray):
        if self.potential_weights_data is not None:
            print(POTENTIAL_WEIGHTS_DATA_ALREADY_SET_WARNING.format(self.atomic_number, self.configuration_id))

        self._check_weights_data(weights_data)
        self.potential_weights_data = weights_data
    

    def set_energy_weights_data(self, weights_data: np.ndarray):
        if self.energy_weights_data is not None:
            print(ENERGY_WEIGHTS_DATA_ALREADY_SET_WARNING.format(self.atomic_number, self.configuration_id))

        self._check_weights_data(weights_data)
        self.energy_weights_data = weights_data


    def _check_weights_data(self, weights_data: np.ndarray):
        assert isinstance(weights_data, np.ndarray), \
            WEIGHTS_DATA_NOT_NUMPY_ARRAY_ERROR.format(weights_data)
        assert weights_data.ndim == 1 or weights_data.ndim == 2, \
            WEIGHTS_DATA_NOT_1D_OR_2D_ARRAY_ERROR.format(weights_data)
        if weights_data.ndim == 1:
            weights_data = weights_data.reshape(-1, 1)
        assert weights_data.shape[0] == self.n_filtered, \
            WEIGHTS_DATA_LENGTH_NOT_EQUAL_TO_N_FILTERED_ERROR.format(weights_data.shape[0], self.n_filtered)


    def clear_weights_data(self):
        self.potential_weights_data = None
        self.energy_weights_data    = None
        

    def _check_inputs_type_value_and_shape(self):
        # Type checks
        assert isinstance(self.atomic_number, int), \
            ATOMIC_NUMBER_NOT_INTEGER_ERROR.format(self.atomic_number)
        assert isinstance(self.n_electrons, float), \
            N_ELECTRONS_NOT_FLOAT_ERROR.format(self.n_electrons)
        assert isinstance(self.configuration_id, int), \
            CONFIGURATION_ID_NOT_INT_ERROR.format(self.configuration_id)
        assert isinstance(self.cutoff_radius, float), \
            CUTOFF_RADIUS_NOT_FLOAT_ERROR.format(self.cutoff_radius)
        assert isinstance(self.cutoff_idx, int), \
            CUTOFF_IDX_NOT_INT_ERROR.format(self.cutoff_idx)
        assert isinstance(self.n_points, int), \
            N_POINTS_NOT_INT_ERROR.format(self.n_points)
        assert isinstance(self.n_filtered, int), \
            N_FILTERED_NOT_INT_ERROR.format(self.n_filtered)
        assert isinstance(self.features_list, list), \
            FEATURES_LIST_NOT_LIST_ERROR.format(self.features_list)

        assert isinstance(self.features_data, np.ndarray), \
            FEATURES_DATA_NOT_NUMPY_ARRAY_ERROR.format(self.features_data)
        assert isinstance(self.quadrature_nodes, np.ndarray), \
            QUADRATURE_NODES_NOT_NUMPY_ARRAY_ERROR.format(self.quadrature_nodes)
        assert isinstance(self.quadrature_nodes_filtered, np.ndarray), \
            QUADRATURE_NODES_FILTERED_NOT_NUMPY_ARRAY_ERROR.format(self.quadrature_nodes_filtered)
        assert isinstance(self.quadrature_weights, np.ndarray), \
            QUADRATURE_WEIGHTS_NOT_NUMPY_ARRAY_ERROR.format(self.quadrature_weights)
        assert isinstance(self.quadrature_weights_filtered, np.ndarray), \
            QUADRATURE_WEIGHTS_NOT_NUMPY_ARRAY_ERROR.format(self.quadrature_weights_filtered)
        assert isinstance(self.scf_xc_data, tuple), \
            SCF_XC_DATA_NOT_TUPLE_ERROR.format(self.scf_xc_data)
        assert isinstance(self.forward_pass_xc_data_list, list), \
            FORWARD_PASS_XC_DATA_LIST_NOT_LIST_ERROR.format(self.forward_pass_xc_data_list)

        # Valud and shape checks
        assert self.features_data.ndim == 2, \
            FEATURES_DATA_NOT_2D_ARRAY_ERROR.format(self.features_data)
        assert self.scf_xc_data[0].ndim == 1, \
            SCF_XC_DATA_NOT_1D_ARRAY_ERROR.format(self.scf_xc_data)
        assert self.scf_xc_data[1].ndim == 1, \
            SCF_XC_DATA_NOT_1D_ARRAY_ERROR.format(self.scf_xc_data)
        assert self.atomic_numbers.ndim == 1, \
            ATOMIC_NUMBERS_NOT_1D_ARRAY_ERROR.format(self.atomic_numbers)
        assert self.quadrature_nodes.ndim == 1, \
            QUADRATURE_NODES_NOT_1D_ARRAY_ERROR.format(self.quadrature_nodes)
        assert self.quadrature_nodes_filtered.ndim == 1, \
            QUADRATURE_NODES_FILTERED_NOT_1D_ARRAY_ERROR.format(self.quadrature_nodes_filtered)
        assert self.quadrature_weights.ndim == 1, \
            QUADRATURE_WEIGHTS_NOT_1D_ARRAY_ERROR.format(self.quadrature_weights)
        assert self.quadrature_weights_filtered.ndim == 1, \
            QUADRATURE_WEIGHTS_NOT_1D_ARRAY_ERROR.format(self.quadrature_weights_filtered)
        assert self.cutoff_radius > 0, \
            CUTOFF_RADIUS_NOT_POSITIVE_ERROR.format(self.cutoff_radius)
        assert self.cutoff_idx >= 0, \
            CUTOFF_IDX_NOT_NON_NEGATIVE_ERROR.format(self.cutoff_idx)
        assert self.n_points > 0, \
            N_POINTS_NOT_POSITIVE_ERROR.format(self.n_points)
        assert self.n_filtered > 0, \
            N_FILTERED_NOT_POSITIVE_ERROR.format(self.n_filtered)
        assert isinstance(self.folder_label, str), \
            FOLDER_LABEL_NOT_STRING_ERROR.format(self.folder_label)
        assert isinstance(self.scf_folder_path, str), \
            SCF_FOLDER_PATH_NOT_STRING_ERROR.format(self.scf_folder_path)
        assert isinstance(self.forward_pass_folder_path_list, list), \
            FORWARD_PASS_FOLDER_PATH_LIST_NOT_LIST_ERROR.format(self.forward_pass_folder_path_list)
        if len(self.forward_pass_xc_data_list) > 0:
            assert self.forward_pass_xc_data_list[0][0].ndim == 1, \
                FORWARD_PASS_XC_DATA_NOT_1D_ARRAY_ERROR.format(self.forward_pass_xc_data_list)
            assert self.forward_pass_xc_data_list[0][1].ndim == 1, \
                FORWARD_PASS_XC_DATA_NOT_1D_ARRAY_ERROR.format(self.forward_pass_xc_data_list)
        for forward_pass_folder_path in self.forward_pass_folder_path_list:
            assert isinstance(forward_pass_folder_path, str), \
                FORWARD_PASS_FOLDER_PATH_NOT_STRING_ERROR.format(forward_pass_folder_path)


    def print_info(self):
        """
        Print information about the single configuration data.
        """
        print(f"{'='*60}")
        print(f"\t\t Single Configuration Data")
        print(f"{'='*60}")

        print(f"Atomic number                      : {self.atomic_number}")
        print(f"Number of electrons                : {self.n_electrons}")
        print(f"Configuration id                   : {self.configuration_id}")
        print(f"Cutoff radius                      : {self.cutoff_radius}")
        print(f"Cutoff index                       : {self.cutoff_idx}")
        print(f"Number of points                   : {self.n_points}")
        print(f"Number of filtered points          : {self.n_filtered}")
        print(f"Folder label                       : {self.folder_label}")
        print(f"Features list                      : {self.features_list}")
        print(f"SCF XC functional                  : {self.scf_xc_functional}")
        print(f"SCF folder path                    : {self.scf_folder_path}")
        print(f"Forward pass XC functional list    : {self.forward_pass_xc_functional_list}")
        print(f"Forward pass folder path list      : {self.forward_pass_folder_path_list}")
        
        print()
        print(f"shape of quadrature_nodes          : Array of shape {self.quadrature_nodes.shape}")
        print(f"shape of quadrature_nodes_filtered : Array of shape {self.quadrature_nodes_filtered.shape}")
        print(f"shape of features_data             : Array of shape {self.features_data.shape}")
        print(f"shape of scf_xc_data               : Tuple of 4 elements")
        print(f"    - v_x: Array of shape {self.scf_xc_data[0].shape}")
        print(f"    - v_c: Array of shape {self.scf_xc_data[1].shape}")
        print( "    - e_x: {}".format(f"Array of shape {self.scf_xc_data[2].shape}" if self.scf_xc_data[2] is not None else "None"))
        print( "    - e_c: {}".format(f"Array of shape {self.scf_xc_data[3].shape}" if self.scf_xc_data[3] is not None else "None"))
        print(f"shape of forward_pass_xc_data_list : List of {len(self.forward_pass_xc_data_list)} tuples, each with 4 elements")
        print(f"    - v_x: Array of shape {self.forward_pass_xc_data_list[0][0].shape}")
        print(f"    - v_c: Array of shape {self.forward_pass_xc_data_list[0][1].shape}")
        print( "    - e_x: {}".format(f"Array of shape {self.forward_pass_xc_data_list[0][2].shape}" if self.forward_pass_xc_data_list[0][2] is not None else "None"))
        print( "    - e_c: {}".format(f"Array of shape {self.forward_pass_xc_data_list[0][3].shape}" if self.forward_pass_xc_data_list[0][3] is not None else "None"))
        print(f"{'='*60}")
        print()


    def get_xc_data(self, functional: str) -> XCDataType:
        if functional == self.scf_xc_functional:
            return self.scf_xc_data
        if functional in self.forward_pass_xc_functional_list:
            return self.forward_pass_xc_data_list[self.forward_pass_xc_functional_list.index(functional)]
        raise ValueError(FUNCTIONAL_NOT_IN_CONFIGURATION_DATA_ERROR.format(functional, self.scf_xc_functional, self.forward_pass_xc_functional_list))


    @staticmethod
    def extract_potential_component(xc_data: XCDataType, target_component: PotentialTargetComponent) -> np.ndarray:
        """
        Extract the component of the XC data (potential: v_x, v_c, v_xc, v_x_v_c).
        """
        v_x, v_c, _, _ = xc_data
        if target_component in ["v_x", "enhancement_factor_x"]:
            return v_x
        elif target_component == "v_c":
            return v_c
        elif target_component == "v_xc":
            return v_x + v_c
        elif target_component == "v_x_v_c":
            return np.column_stack([v_x, v_c])
        else:
            raise ValueError(TARGET_COMPONENT_NOT_IN_VALID_LIST_ERROR.format(VALID_POTENTIAL_TARGET_COMPONENTS, target_component))


    @staticmethod
    def extract_energy_component(xc_data: XCDataType, target_component: EnergyTargetComponent) -> np.ndarray:
        """
        Extract the energy density component from XC data.
        target_component: "e_xc", "e_x", "e_c", or "e_x_e_c".
        Returns 1d array (e_xc/e_x/e_c) or 2d array (e_x_e_c).
        """
        _, _, e_x, e_c = xc_data
        if e_x is None or e_c is None:
            raise ValueError("XC data does not include energy density (e_x, e_c).")
        if target_component in ["e_x", "enhancement_factor_x"]:
            return e_x
        elif target_component == "e_c":
            return e_c
        elif target_component == "e_xc":
            return e_x + e_c
        elif target_component == "e_x_e_c":
            return np.column_stack([e_x, e_c])
        else:
            raise ValueError(TARGET_COMPONENT_NOT_IN_VALID_LIST_ERROR.format(VALID_ENERGY_TARGET_COMPONENTS, target_component))



    def extract_potential_target(
        self, 
        target_component     : PotentialTargetComponent,
        target_functional    : str,
        reference_functional : Optional[str] = None,
    ) -> np.ndarray:
        """
        Extract the potential target from the XC data.
        """
        v_xc_target = self.extract_potential_component(self.get_xc_data(target_functional), target_component)
        if reference_functional is not None:
            v_xc_reference = self.extract_potential_component(self.get_xc_data(reference_functional), target_component)
            potential_target = v_xc_target - v_xc_reference
        else:
            potential_target = v_xc_target

        if target_component == "enhancement_factor_x":
            from ..xc.lda import lda_exchange_potential_generic
            lda_exchange_potential = lda_exchange_potential_generic(self.rho)
            # Safeguard: at large radius, v_lda → 0 causes division by zero, Inf, or NaN.
            # When |v_lda| is too small, F_x = v_x/v_lda is numerically unreliable (0/0 or underflow).
            # Use F_x = 1.0 (LDA limit) as fallback when denominator is tiny.
            mask_safe = np.abs(lda_exchange_potential) > 1e-15
            potential_target = np.where(
                mask_safe,
                potential_target / lda_exchange_potential,
                1.0,  # F_x → 1 at ρ→0 (LDA limit)
            )
        else:
            pass
        
        return potential_target


    def extract_energy_target(
        self, 
        target_component     : EnergyTargetComponent,
        target_functional    : str,
        reference_functional : Optional[str] = None,
    ) -> np.ndarray:
        """
        Extract the energy target from the XC data.
        """
        e_xc_target = self.extract_energy_component(self.get_xc_data(target_functional), target_component)
        
        # If reference functional is provided, compute the delta energy target
        if reference_functional is not None:
            e_xc_reference = self.extract_energy_component(self.get_xc_data(reference_functional), target_component)
            energy_target = e_xc_target - e_xc_reference
        else:
            energy_target = e_xc_target
        
        if target_component == "enhancement_factor_x":
            from ..xc.lda import lda_exchange_energy_density_generic
            lda_exchange_energy_density = lda_exchange_energy_density_generic(self.rho) * self.rho
            # Safeguard: at large radius, ε_lda*ρ → 0 causes division by zero or numerical instability.
            # When denominator is too small, F_x = e_x/(ε_lda*ρ) is unreliable (0/0 or underflow → 0).
            # Use F_x = 1.0 (LDA limit) as fallback when denominator is tiny.
            mask_safe = np.abs(lda_exchange_energy_density) > 1e-20
            energy_target = np.where(
                mask_safe,
                energy_target / lda_exchange_energy_density,
                1.0,  # F_x → 1 at ρ→0 (LDA limit)
            )
        else:
            pass
        
        return energy_target


    def compute_total_energy_true(
        self,
        target_component     : EnergyTargetComponent,
        target_functional    : str,
        reference_functional : Optional[str] = None,
    ) -> Union[float, np.ndarray]:
        """
        Compute total energy E_true = ∫ 4πr² w·e for total_energy loss.
        Uses raw (untransformed) energy density from XC data.
        All quadrature uses filtered nodes/weights to avoid truncation error.
        """
        e_target = self.extract_energy_component(self.get_xc_data(target_functional), target_component)
        if reference_functional is not None:
            e_reference = self.extract_energy_component(self.get_xc_data(reference_functional), target_component)
            energy_density = e_target - e_reference
        else:
            energy_density = e_target

        r = np.asarray(self.quadrature_nodes_filtered).reshape(-1, 1)
        w = np.asarray(self.quadrature_weights_filtered).reshape(-1, 1)
        e = np.asarray(energy_density, dtype=float)
        if e.ndim == 1:
            e = e.reshape(-1, 1)
        integrand = 4 * np.pi * (r ** 2) * w * e
        result = np.sum(integrand, axis=0)
        return result[0] if result.size == 1 else result


    def get_features(self, feature_name: str) -> np.ndarray:
        return self.features_data[:, self.features_list.index(feature_name)]


    @property
    def atomic_numbers(self) -> np.ndarray:
        return np.full(self.n_filtered, self.atomic_number)


    @property
    def rho(self) -> np.ndarray:
        return self.features_data[:, self.features_list.index("rho")]


    @property
    def rho_filtered(self) -> np.ndarray:
        return self.features_data[:, self.features_list.index("rho")][:self.n_filtered]


    @property
    def quadrature_nodes_filtered(self) -> np.ndarray:
        return self.quadrature_nodes[:self.n_filtered]

    @property
    def quadrature_weights_filtered(self) -> np.ndarray:
        return self.quadrature_weights[:self.n_filtered]



class DataLoader:

    @staticmethod
    def check_and_normalize_features_list(features_list: List[str]) -> List[str]:
        """
        Normalize feature names to the valid features list and validate them.
        """
        if not isinstance(features_list, list):
            raise ValueError(FEATURES_LIST_NOT_LIST_ERROR.format(type(features_list)))
        normalized = []
        for feature in features_list:
            if not isinstance(feature, str):
                raise ValueError(
                    format_invalid_feature_error(
                        feature,
                        NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL,
                    )
                )
            normalized_feature = FEATURE_ALIASES.get(feature, feature)
            if normalized_feature not in NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL:
                raise ValueError(
                    format_invalid_feature_error(
                        feature,
                        NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL,
                    )
                )
            normalized.append(normalized_feature)
        return normalized


    @staticmethod
    def _rho_is_monotone_nonincreasing(rho: np.ndarray, rtol: float = 1e-10) -> bool:
        """
        True if rho[i] >= rho[i+1] (up to relative tolerance) along the grid; expected for smoothed
        spherical density when radial nodes increase outward.
        """
        r = np.asarray(rho, dtype=float).reshape(-1)
        if r.size < 2:
            return True
        step = r[:-1] - r[1:]
        scale = np.maximum(np.abs(r[:-1]), np.abs(r[1:]))
        scale = np.maximum(scale, 1e-300)
        return bool(np.all(step >= -rtol * scale))


    @staticmethod
    def load_quadrature_nodes_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load quadrature nodes data from a folder.
        
        Args:
            folder_path            : Path to folder containing quadrature nodes data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : Quadrature nodes data
        """

        quadrature_nodes = np.loadtxt(os.path.join(folder_path, "quadrature_nodes.txt"))

        # Basic type and shape checks
        if not isinstance(quadrature_nodes, np.ndarray):
            raise ValueError(QUADRATURE_NODES_NOT_NUMPY_ARRAY_ERROR.format(quadrature_nodes))
        if not quadrature_nodes.ndim == 1:
            raise ValueError(QUADRATURE_NODES_NOT_1D_ARRAY_ERROR.format(quadrature_nodes))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(quadrature_nodes) != reference_array_length:
                raise ValueError(QUADRATURE_NODES_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(quadrature_nodes), reference_array_length))
            
        return quadrature_nodes


    @staticmethod
    def load_quadrature_weights_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load quadrature weights data from a folder.
        
        Args:
            folder_path            : Path to folder containing quadrature weights data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : Quadrature weights data
        """

        quadrature_weights = np.loadtxt(os.path.join(folder_path, "quadrature_weights.txt"))

        # Basic type and shape checks
        if not isinstance(quadrature_weights, np.ndarray):
            raise ValueError(QUADRATURE_WEIGHTS_NOT_NUMPY_ARRAY_ERROR.format(quadrature_weights))
        if not quadrature_weights.ndim == 1:
            raise ValueError(QUADRATURE_WEIGHTS_NOT_1D_ARRAY_ERROR.format(quadrature_weights))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(quadrature_weights) != reference_array_length:
                raise ValueError(QUADRATURE_WEIGHTS_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(quadrature_weights), reference_array_length))
            
        return quadrature_weights


    @staticmethod
    def load_rho_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load rho data from a folder.
        
        Args:
            folder_path            : Path to folder containing rho data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : Rho data
        """
        rho = np.loadtxt(os.path.join(folder_path, "rho.txt"))
        
        # Basic type and shape checks
        if not isinstance(rho, np.ndarray):
            raise ValueError(RHO_NOT_NUMPY_ARRAY_ERROR.format(rho))
        if not rho.ndim == 1:
            raise ValueError(RHO_NOT_1D_ARRAY_ERROR.format(rho))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(rho) != reference_array_length:
                raise ValueError(RHO_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(rho), reference_array_length))
            
        return rho


    @staticmethod
    def load_grad_rho_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
        use_feature_round_off  : bool = True,
    ) -> np.ndarray:
        """
        Load grad rho data from a folder.
        
        Args:
            folder_path            : Path to folder containing grad rho data    
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : Grad rho data
        """
        grad_rho = np.loadtxt(os.path.join(folder_path, "grad_rho.txt"))
        
        # Basic type and shape checks
        if not isinstance(grad_rho, np.ndarray):
            raise ValueError(GRAD_RHO_NOT_NUMPY_ARRAY_ERROR.format(grad_rho))
        if not grad_rho.ndim == 1:
            raise ValueError(GRAD_RHO_NOT_1D_ARRAY_ERROR.format(grad_rho))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(grad_rho) != reference_array_length:
                raise ValueError(GRAD_RHO_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(grad_rho), reference_array_length))
            
        if use_feature_round_off:
            # Round off to the upper bound
            grad_rho = np.minimum(grad_rho, GRAD_RHO_UPPER_BOUND)

        return grad_rho


    @staticmethod
    def load_lap_rho_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load lap rho data from a folder.
        
        Args:
            folder_path            : Path to folder containing lap rho data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : Lap rho data
        """
        lap_rho = np.loadtxt(os.path.join(folder_path, "lap_rho.txt"))
        
        # Basic type and shape checks
        if not isinstance(lap_rho, np.ndarray):
            raise ValueError(LAP_RHO_NOT_NUMPY_ARRAY_ERROR.format(lap_rho))
        if not lap_rho.ndim == 1:
            raise ValueError(LAP_RHO_NOT_1D_ARRAY_ERROR.format(lap_rho))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(lap_rho) != reference_array_length:
                raise ValueError(LAP_RHO_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(lap_rho), reference_array_length))
            
        # No round off needed for lap_rho
        return lap_rho


    @staticmethod
    def load_hartree_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
        use_feature_round_off  : bool = True,
    ) -> np.ndarray:
        """
        Load hartree data from a folder.
        
        Args:
            folder_path            : Path to folder containing hartree data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : Hartree data
        """
        hartree = np.loadtxt(os.path.join(folder_path, "hartree.txt"))
        
        # Basic type and shape checks
        if not isinstance(hartree, np.ndarray):
            raise ValueError(HARTREE_NOT_NUMPY_ARRAY_ERROR.format(hartree))
        if not hartree.ndim == 1:
            raise ValueError(HARTREE_NOT_1D_ARRAY_ERROR.format(hartree))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(hartree) != reference_array_length:
                raise ValueError(HARTREE_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(hartree), reference_array_length))
            
        if use_feature_round_off:
            # Round off to the lower bound
            hartree = np.maximum(hartree, HARTREE_LOWER_BOUND)
        
        return hartree


    @staticmethod
    def compute_lda_xc_potential_data(
        rho : np.ndarray,
        use_feature_round_off : bool = True,
    ) -> np.ndarray:
        """
        Compute the XC potential of LDA functional for given rho data.
        
        Args:
            rho : Rho data
        
        Returns:
            np.ndarray of shape (n_points,) : Exchange and correlation potential data
        """
        # Check if rho is a numpy array
        assert isinstance(rho, np.ndarray), RHO_NOT_NUMPY_ARRAY_ERROR.format(rho)
        assert rho.ndim == 1, RHO_NOT_1D_ARRAY_ERROR.format(rho)
        
        # Compute LDA exchange and correlation
        lda_exchange = lda_exchange_potential_generic(rho)
        lda_correlation = lda_correlation_potential_generic(rho)
        lda_xc_potential = lda_exchange + lda_correlation
    
        if use_feature_round_off:
            # Round off to the upper bound
            lda_xc_potential = np.minimum(lda_xc_potential, LDA_XC_POTENTIAL_UPPER_BOUND)
        
        return lda_xc_potential

    @staticmethod
    def compute_grad_rho_norm(grad_rho: np.ndarray) -> np.ndarray:
        return np.abs(grad_rho)

    @staticmethod
    def compute_grad_rho_reduced(
        rho           : np.ndarray,
        grad_rho      : np.ndarray,
        rho_threshold : float = 0.0,
    ) -> np.ndarray:
        """s = |∇ρ| / (2 * kf * ρ^(4/3)). rho_threshold: lower bound for rho to avoid division issues."""
        rho_safe = np.maximum(rho, rho_threshold)
        kf = (3 * np.pi**2) ** (1 / 3)
        return np.abs(grad_rho) / (2 * kf * (rho_safe ** (4 / 3)))

    @staticmethod
    def compute_lap_rho_reduced(rho: np.ndarray, lap_rho: np.ndarray) -> np.ndarray:
        kf = (3 * np.pi**2) ** (1 / 3)
        return lap_rho / (4 * (kf ** 2) * (rho ** (5 / 3)))
    

    @staticmethod
    def _atomic_solver_from_output_in_folder(folder_path: str) -> "AtomicDFTSolver":
        from ..solver import AtomicDFTSolver

        out_file_path = os.path.join(folder_path, "out.txt")
        if not os.path.exists(out_file_path):
            parent_path = os.path.dirname(folder_path)
            out_file_path = os.path.join(parent_path, "out.txt")
        if not os.path.exists(out_file_path):
            raise FileNotFoundError(
                f"out.txt not found in {folder_path} or {os.path.dirname(folder_path)}, "
                "required for solver-based density derivatives"
            )
        return AtomicDFTSolver.from_output_file(out_file_path, verbose=False)

    @classmethod
    def load_features(
        cls,
        folder_path                : str,
        features_list              : List[str],
        reference_array_length     : Optional[int] = None,  # For debugging and error handling
        use_feature_round_off      : bool = True,
        smooth_reduced_laplacian   : bool = False,
        smooth_rho_for_derivatives : bool = False,
        smooth_rho_r_a             : float = 7.0,
        smooth_rho_r_b             : float = 10.0,
        smooth_rho_r_c             : float = 13.0,
        smooth_rho_log_floor       : float = 1e-30,
    ) -> np.ndarray:
        """
        Load features from a folder.
        
        Args:
            folder_path              : Path to folder containing features
            features_list            : List of features to load
            reference_array_length   : Length of the reference array to check against. For debugging and error handling.
            use_feature_round_off    : Whether to round off the features to the upper or lower bound
            smooth_reduced_laplacian : Whether to smooth the reduced Laplacian at large radius
                (uses out.txt in folder_path to construct solver for weak-form Laplacian)
            smooth_rho_for_derivatives : If True, apply log-linear asymptotic smoothing to rho, then compute
                grad_rho / lap_rho (and reduced s, q) from that rho with the FEM operators in out.txt instead
                of reading grad_rho.txt / lap_rho.txt.
            smooth_rho_r_a, smooth_rho_r_b, smooth_rho_r_c :
                Radii in Bohr with ``r_a < r_b < r_c``; see ``DataProcessor.smooth_rho_with_loglinear_asymptote``.
            smooth_rho_log_floor : Floor for log fit inside the asymptotic window.

        Returns:
            np.ndarray of shape (n_points, len(features_list)) : Features data
        """
        # Normalize features list, and rho data must be included
        features_list = cls.check_and_normalize_features_list(features_list)
        if "rho" not in features_list:
            raise ValueError(RHO_DATA_MUST_BE_INCLUDED_IN_FEATURES_LIST_ERROR)

        feature_data_list : List[np.ndarray] = []
        rho = cls.load_rho_data(folder_path, reference_array_length)

        if smooth_rho_for_derivatives:
            r_grid = cls.load_quadrature_nodes_data(folder_path, reference_array_length)
            if len(r_grid) != len(rho):
                raise ValueError(
                    f"Quadrature nodes length {len(r_grid)} does not match rho length {len(rho)} "
                    f"in {folder_path} (required for asymptotic rho smoothing)."
                )
            rho = DataProcessor.smooth_rho_with_loglinear_asymptote(
                r_grid,
                rho,
                smooth_rho_r_a,
                smooth_rho_r_b,
                smooth_rho_r_c,
                rho_floor=smooth_rho_log_floor,
            )

        lap_rho_strong = None  # matches lap_rho.txt convention (FEM strong form)
        lap_rho_weak = None    # weak-form lap for lap_rho_reduced when smooth_reduced_laplacian
        rho_for_features = np.maximum(rho, RHO_LOWER_BOUND) if use_feature_round_off else rho
        grad_rho = None
        solver = None
        density_calculator = None

        def ensure_solver() -> None:
            nonlocal solver, density_calculator
            if solver is None:
                solver = cls._atomic_solver_from_output_in_folder(folder_path)
                density_calculator = solver.density_calculator

        for feature in features_list:
            if feature == "rho":
                feature_data_list.append(rho_for_features)
            elif feature == "grad_rho":
                if smooth_rho_for_derivatives:
                    ensure_solver()
                    if grad_rho is None:
                        grad_rho = density_calculator.compute_density_gradient(rho)
                        if use_feature_round_off:
                            grad_rho = np.minimum(grad_rho, GRAD_RHO_UPPER_BOUND)
                else:
                    if grad_rho is None:
                        grad_rho = cls.load_grad_rho_data(
                            folder_path,
                            reference_array_length,
                            use_feature_round_off=use_feature_round_off,
                        )
                feature_data_list.append(grad_rho)
            elif feature == "grad_rho_norm":
                if smooth_rho_for_derivatives:
                    ensure_solver()
                    if grad_rho is None:
                        grad_rho = density_calculator.compute_density_gradient(rho)
                        if use_feature_round_off:
                            grad_rho = np.minimum(grad_rho, GRAD_RHO_UPPER_BOUND)
                else:
                    if grad_rho is None:
                        grad_rho = cls.load_grad_rho_data(
                            folder_path,
                            reference_array_length,
                            use_feature_round_off=use_feature_round_off,
                        )
                feature_data_list.append(cls.compute_grad_rho_norm(grad_rho))
            elif feature == "lap_rho":
                if smooth_rho_for_derivatives:
                    ensure_solver()
                    if lap_rho_strong is None:
                        lap_rho_strong = density_calculator.compute_density_laplacian(rho)
                else:
                    if lap_rho_strong is None:
                        lap_rho_strong = cls.load_lap_rho_data(folder_path, reference_array_length)
                feature_data_list.append(lap_rho_strong)
            elif feature == "hartree":
                feature_data_list.append(
                    cls.load_hartree_data(
                        folder_path,
                        reference_array_length,
                        use_feature_round_off=use_feature_round_off,
                    )
                )
            elif feature == "lda_xc_potential":
                feature_data_list.append(
                    cls.compute_lda_xc_potential_data(
                        rho_for_features,
                        use_feature_round_off=use_feature_round_off,
                    )
                )
            elif feature == "grad_rho_reduced":
                if smooth_rho_for_derivatives:
                    ensure_solver()
                    if grad_rho is None:
                        grad_rho = density_calculator.compute_density_gradient(rho)
                        if use_feature_round_off:
                            grad_rho = np.minimum(grad_rho, GRAD_RHO_UPPER_BOUND)
                else:
                    if grad_rho is None:
                        grad_rho = cls.load_grad_rho_data(
                            folder_path,
                            reference_array_length,
                            use_feature_round_off=use_feature_round_off,
                        )
                feature_data_list.append(cls.compute_grad_rho_reduced(rho, grad_rho))
            elif feature == "lap_rho_reduced":
                from .data_generation import DataGenerator
                if smooth_rho_for_derivatives:
                    ensure_solver()
                    if smooth_reduced_laplacian:
                        if lap_rho_weak is None:
                            lap_rho_weak = DataGenerator.calculate_laplacian_in_weak_form(rho, solver)
                        lap_for_reduced = lap_rho_weak
                    else:
                        if lap_rho_strong is None:
                            lap_rho_strong = density_calculator.compute_density_laplacian(rho)
                        lap_for_reduced = lap_rho_strong
                else:
                    if smooth_reduced_laplacian:
                        wk_solver = cls._atomic_solver_from_output_in_folder(folder_path)
                        lap_for_reduced = DataGenerator.calculate_laplacian_in_weak_form(rho, wk_solver)
                    else:
                        if lap_rho_strong is None:
                            lap_rho_strong = cls.load_lap_rho_data(folder_path, reference_array_length)
                        lap_for_reduced = lap_rho_strong
                feature_data_list.append(cls.compute_lap_rho_reduced(rho, lap_for_reduced))
            else:
                raise ValueError(format_invalid_feature_error(feature, NORMALIZED_VALID_FEATURES_LIST_FOR_POTENTIAL))
        
        # Stack features data: column_stack creates (n_points, n_features) shape
        features_data = np.column_stack(feature_data_list)
        
        # Shape checks: features_data.shape should be (n_points, n_features)
        # So features_data.shape[0] is n_points, features_data.shape[1] is n_features (len(features_list))
        if features_data.shape[1] != len(features_list):
            raise ValueError(FEATURES_DATA_LENGTH_NOT_EQUAL_TO_FEATURES_LIST_LENGTH_ERROR.format(features_data.shape[1], len(features_list)))
        if reference_array_length is not None:
            if features_data.shape[0] != reference_array_length:
                raise ValueError(FEATURES_DATA_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(features_data.shape[0], reference_array_length))
        return features_data


    @staticmethod
    def load_v_x_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load v_x_pbe0 data from a folder.
        
        Args:
            folder_path            : Path to folder containing v_x data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : V_x data
        """
        # Try to load v_x.txt first
        v_x_file = os.path.join(folder_path, "v_x.txt")
        if os.path.exists(v_x_file):
            file_path = v_x_file
        else:
            # Look for alternative files like v_x_pbe0.txt, but exclude *_uniform.txt
            pattern = os.path.join(folder_path, "v_x_*.txt")
            matching_files = [f for f in glob.glob(pattern) if not f.endswith("_uniform.txt")]
            if not matching_files:
                raise FileNotFoundError(V_X_FILE_NOT_FOUND_ERROR.format(folder_path))
            # Use the first matching file (sorted for consistency)
            file_path = sorted(matching_files)[0]
        
        v_x = np.loadtxt(file_path)
        
        # Basic type and shape checks
        if not isinstance(v_x, np.ndarray):
            raise ValueError(V_X_NOT_NUMPY_ARRAY_ERROR.format(v_x))
        if not v_x.ndim == 1:
            raise ValueError(V_X_NOT_1D_ARRAY_ERROR.format(v_x))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(v_x) != reference_array_length:
                raise ValueError(V_X_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(v_x), reference_array_length))
            
        return v_x


    @staticmethod
    def load_v_c_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load v_c data from a folder.
        
        Args:
            folder_path            : Path to folder containing v_c data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : V_c data
        """
        # Try to load v_c.txt first
        v_c_file = os.path.join(folder_path, "v_c.txt")
        if os.path.exists(v_c_file):
            file_path = v_c_file
        else:
            # Look for alternative files like v_c_pbe0.txt, but exclude *_uniform.txt
            pattern = os.path.join(folder_path, "v_c_*.txt")
            matching_files = [f for f in glob.glob(pattern) if not f.endswith("_uniform.txt")]
            if not matching_files:
                raise FileNotFoundError(V_C_FILE_NOT_FOUND_ERROR.format(folder_path))
            # Use the first matching file (sorted for consistency)
            file_path = sorted(matching_files)[0]
        
        v_c = np.loadtxt(file_path)
        
        # Basic type and shape checks
        if not isinstance(v_c, np.ndarray):
            raise ValueError(V_C_NOT_NUMPY_ARRAY_ERROR.format(v_c))
        if not v_c.ndim == 1:
            raise ValueError(V_C_NOT_1D_ARRAY_ERROR.format(v_c))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(v_c) != reference_array_length:
                raise ValueError(V_C_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(v_c), reference_array_length))
            
        return v_c


    @staticmethod
    def load_e_x_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load e_x data from a folder.
        
        Args:
            folder_path            : Path to folder containing e_x data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : E_x data
        """
        # Try to load e_x.txt first
        e_x_file = os.path.join(folder_path, "e_x.txt")
        if os.path.exists(e_x_file):
            file_path = e_x_file
        else:
            # Look for alternative files like e_x_pbe0.txt, but exclude *_uniform.txt
            pattern = os.path.join(folder_path, "e_x_*.txt")
            matching_files = [f for f in glob.glob(pattern) if not f.endswith("_uniform.txt")]
            if not matching_files:
                raise FileNotFoundError(E_X_FILE_NOT_FOUND_ERROR.format(folder_path))
            # Use the first matching file (sorted for consistency)
            file_path = sorted(matching_files)[0]
        
        e_x = np.loadtxt(file_path)
        
        # Basic type and shape checks
        if not isinstance(e_x, np.ndarray):
            raise ValueError(E_X_NOT_NUMPY_ARRAY_ERROR.format(e_x))
        if not e_x.ndim == 1:
            raise ValueError(E_X_NOT_1D_ARRAY_ERROR.format(e_x))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(e_x) != reference_array_length:
                raise ValueError(E_X_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(e_x), reference_array_length))
            
        return e_x


    @staticmethod
    def load_e_c_data(
        folder_path            : str,
        reference_array_length : Optional[int] = None,  # For debugging and error handling
    ) -> np.ndarray:
        """
        Load e_c data from a folder.
        

        Args:
            folder_path            : Path to folder containing e_c data
            reference_array_length : Length of the reference array to check against. For debugging and error handling.
        
        Returns:
            np.ndarray of shape (n_points,) : E_c data
        """
        # Try to load e_c.txt first
        e_c_file = os.path.join(folder_path, "e_c.txt")
        if os.path.exists(e_c_file):
            file_path = e_c_file
        else:
            # Look for alternative files like e_c_pbe0.txt, but exclude *_uniform.txt
            pattern = os.path.join(folder_path, "e_c_*.txt")
            matching_files = [f for f in glob.glob(pattern) if not f.endswith("_uniform.txt")]
            if not matching_files:
                raise FileNotFoundError(E_C_FILE_NOT_FOUND_ERROR.format(folder_path))
            # Use the first matching file (sorted for consistency)
            file_path = sorted(matching_files)[0]
        
        e_c = np.loadtxt(file_path)
        
        # Basic type and shape checks
        if not isinstance(e_c, np.ndarray):
            raise ValueError(E_C_NOT_NUMPY_ARRAY_ERROR.format(e_c))
        if not e_c.ndim == 1:
            raise ValueError(E_C_NOT_1D_ARRAY_ERROR.format(e_c))
        
        # Extra checks for reference array length
        if reference_array_length is not None:
            if not isinstance(reference_array_length, int):
                raise ValueError(REFERENCE_ARRAY_LENGTH_NOT_INTEGER_ERROR.format(reference_array_length))
            if len(e_c) != reference_array_length:
                raise ValueError(E_C_LENGTH_NOT_EQUAL_TO_REFERENCE_ARRAY_LENGTH_ERROR.format(len(e_c), reference_array_length))
            
        return e_c


    @staticmethod
    def _find_configuration_folder(
        data_root  : str,
        config_idx : int,
    ) -> Optional[str]:
        """
        Find configuration folder with backward compatibility.
        
        Tries to find folder in this order:
        1. configuration_XXX (new format)
        2. atom_XXX (old format, for backward compatibility)
        
        Args:
            data_root  : Root directory of the dataset
            config_idx : Configuration index
        
        Returns:
            str: Path to configuration folder, or None if not found
        """
        # Try new format first
        config_folder = os.path.join(data_root, f"configuration_{config_idx:03d}")
        if os.path.exists(config_folder) and os.path.isdir(config_folder):
            return config_folder
        
        # Fallback to old format for backward compatibility
        atom_folder = os.path.join(data_root, f"atom_{config_idx:03d}")
        if os.path.exists(atom_folder) and os.path.isdir(atom_folder):
            return atom_folder
        
        return None

    @staticmethod
    def _load_metadata_from_meta(
        config_folder : str,
        config_idx    : int,
    ) -> Tuple[int, Optional[float], Optional[bool]]:
        """
        Load atomic_number, n_electrons, and converged status from meta.json file.
        
        If meta.json is not found or cannot be read, fallback to backward compatibility mode:
        - Check if folder name starts with 'atom_'
        - Extract atomic_number from folder name (atom_XXX -> XXX)
        - Assume n_electrons = atomic_number (neutral atom)
        - Assume converged = True
        
        Args:
            config_folder : Path to configuration folder (e.g., data_root/configuration_001 or data_root/atom_001)
            config_idx    : Configuration index (extracted from folder name)
        
        Returns:
            Tuple[int, Optional[float], Optional[bool]]: (atomic_number, n_electrons, converged)
                - atomic_number: from meta.json, or extracted from folder name (backward compatibility)
                - n_electrons: from meta.json, or atomic_number (backward compatibility for neutral atoms)
                - converged: from meta.json, or True (backward compatibility)
        
        Raises:
            UserWarning: If meta.json is not found or cannot be read, uses backward compatibility mode
        """
        meta_path = os.path.join(config_folder, "meta.json")
        
        if os.path.exists(meta_path):
            try:
                with open(meta_path, "r") as meta_file:
                    meta_data = json.load(meta_file)
                
                # Load atomic_number (with fallback to config_idx)
                atomic_number = int(meta_data.get("atomic_number", 0))
                if atomic_number <= 0:
                    atomic_number = config_idx
                    print(ATOMIC_NUMBER_NOT_FOUND_IN_META_WARNING.format(config_folder, config_idx))
                
                # Load n_electrons
                n_electrons = meta_data.get("n_electrons")
                if n_electrons is not None:
                    n_electrons = float(n_electrons)
                
                # Load converged
                converged = meta_data.get("converged")
                if converged is not None:
                    converged = bool(converged)
                
                return atomic_number, n_electrons, converged
            except Exception:
                pass
        
        # Backward compatibility: if meta.json doesn't exist, try to extract from folder name
        folder_name = os.path.basename(config_folder)
        
        # Check if folder name starts with 'atom_' (old format)
        if folder_name.startswith("atom_"):
            try:
                # Extract atomic_number from folder name (atom_XXX -> XXX)
                atomic_number_str = folder_name.replace("atom_", "")
                atomic_number = int(atomic_number_str)
                
                # For backward compatibility: assume neutral atom (n_electrons = atomic_number)
                n_electrons = float(atomic_number)
                
                # For backward compatibility: assume converged
                converged = True
                
                return atomic_number, n_electrons, converged
            except (ValueError, AttributeError):
                # If extraction fails, folder name format is invalid (e.g., atom_abc instead of atom_010)
                print(ATOMIC_NUMBER_EXTRACTION_FAILED_WARNING.format(folder_name, config_idx))
                return config_idx, None, None
        
        # If folder name doesn't start with 'atom_', it's likely a new format (configuration_XXX)
        # but meta.json is missing, so we can't extract atomic_number
        print(META_JSON_NOT_FOUND_WARNING.format(folder_name, config_idx))
        return config_idx, None, None

    @staticmethod
    def load_derivative_matrix(
        config_folder_path : str,
        data_root          : str,
    ) -> Tuple[bool, Optional[np.ndarray]]:
        """
        Load derivative matrix information from meta.json and return whether it's shared and the matrix itself.
        
        Args:
            config_folder_path : Path to configuration folder (e.g., data_root/configuration_001 or data_root/atom_001)
            data_root          : Root directory of the dataset
        
        Returns:
            Tuple[bool, Optional[np.ndarray]]:
                - use_shared: True if derivative matrix is shared (stored at dataset root), False if unique per configuration
                - derivative_matrix: None if use_shared=True, otherwise the unique derivative matrix for this configuration
        """
        meta_path = os.path.join(config_folder_path, "meta.json")
        
        if not os.path.exists(meta_path):
            return False, None
        
        try:
            with open(meta_path, "r") as meta_file:
                meta_data = json.load(meta_file)
            
            if "derivative_matrix" not in meta_data:
                return False, None
            
            dm_info = meta_data["derivative_matrix"]
            use_shared = dm_info.get("use_shared", False)
            
            if use_shared:
                # Shared derivative matrix is stored at dataset root, return None here
                # The actual matrix will be loaded at dataset level
                return True, None
            else:
                # Load unique derivative matrix for this atom
                unique_path = os.path.join(data_root, dm_info.get("path", ""))
                if os.path.exists(unique_path):
                    derivative_matrix = np.load(unique_path)
                    return False, derivative_matrix
                else:
                    return False, None
        except Exception:
            # If loading fails, return False, None
            return False, None



    @staticmethod
    def compute_cutoff_index(
        rho           : np.ndarray,
        rho_threshold : float = 1e-6,
        v_x           : Optional[np.ndarray] = None,
        v_c           : Optional[np.ndarray] = None,
        v_x_threshold : Optional[float] = None,
        v_c_threshold : Optional[float] = None,
    ) -> int:
        """
        Compute the cutoff data for given rho, v_x, and v_c data.
        
        Args:
            rho           : Rho data
            rho_threshold : Threshold for rho data
            v_x           : Exchange potential data
            v_c           : Correlation potential data
            v_x_threshold : Threshold for exchange potential data
            v_c_threshold : Threshold for correlation potential data

        Returns:
            int : Cutoff index
        """
        # Check if rho is a numpy array
        assert isinstance(rho, np.ndarray), RHO_NOT_NUMPY_ARRAY_ERROR.format(rho)
        assert rho.ndim == 1, RHO_NOT_1D_ARRAY_ERROR.format(rho)
        
        # Check if v_x is a numpy array
        if v_x is not None:
            assert isinstance(v_x, np.ndarray), V_X_NOT_NUMPY_ARRAY_ERROR.format(v_x)
            assert v_x.ndim == 1, V_X_NOT_1D_ARRAY_ERROR.format(v_x)
        
        # Check if v_c is a numpy array
        if v_c is not None:
            assert isinstance(v_c, np.ndarray), V_C_NOT_NUMPY_ARRAY_ERROR.format(v_c)
            assert v_c.ndim == 1, V_C_NOT_1D_ARRAY_ERROR.format(v_c)
        
        # Create mask based on filtering criteria
        mask = (np.abs(rho) > rho_threshold)
        if v_x is not None:
            mask &= (np.abs(v_x) > v_x_threshold)
        if v_c is not None:
            mask &= (np.abs(v_c) > v_c_threshold)
        
        # Find the cutoff index: the last index in a continuous sequence starting from 0
        # Check if index 0 is valid
        if not mask[0]:
            raise ValueError(CUTOFF_INDEX_NOT_FOUND_ERROR)
        
        # Find the last index in a continuous sequence starting from 0
        cutoff_idx = 0
        for i in range(len(mask)):
            if mask[i]:
                cutoff_idx = i + 1
            else:
                break
        
        return cutoff_idx


    @classmethod
    def load_single_configuration_data(
        cls,
        # Basic arguments
        atomic_number                   : int,        # Required only for documentation and error handling
        n_electrons                     : float,
        features_list                   : List[str],
        scf_xc_functional               : str,
        forward_pass_xc_functional_list : List[str],

        # path control
        scf_folder_path                 : str,
        forward_pass_folder_path_list   : List[str],

        # Extra settings after loading data
        use_radius_cutoff               : bool,
        use_feature_round_off           : bool,
        include_energy_density          : bool,
        folder_label                    : str,

        # Super-parameters for cutoff filtering, valid only when use_radius_cutoff is True
        radius_cutoff_rho_threshold     : Optional[float] = None,
        radius_cutoff_v_x_threshold     : Optional[float] = None,
        radius_cutoff_v_c_threshold     : Optional[float] = None,
        smooth_reduced_laplacian        : bool = False,
        smooth_rho_for_derivatives      : bool = False,
        smooth_rho_r_a                  : float = 7.0,
        smooth_rho_r_b                  : float = 10.0,
        smooth_rho_r_c                  : float = 13.0,
        smooth_rho_log_floor            : float = 1e-30,
    ) -> Optional[SingleConfigurationData]:
        """
        Helper function to load data from a single folder (main or intermediate).
        
        Args:
            atomic_number                 : Atomic number
            n_electrons                   : Number of electrons
            features_list                 : list of features to load
            scf_folder_path               : Path to SCF folder
            forward_pass_folder_path_list : List of path(s) to forward pass folders
            use_radius_cutoff             : Whether to apply cutoff filtering
            use_feature_round_off         : Whether to apply feature round-off (lower/upper bounds)
            folder_label                  : Label for this folder (e.g., "main", "outer_iter_01")
            radius_cutoff_rho_threshold   : Threshold for rho data, valid only when use_radius_cutoff is True
            radius_cutoff_v_x_threshold   : Threshold for v_x data, valid only when use_radius_cutoff is True
            radius_cutoff_v_c_threshold   : Threshold for v_c data, valid only when use_radius_cutoff is True

        Returns:
            SingleConfigurationData object
        """
        try:
            # Check if folders exist
            if not os.path.exists(scf_folder_path):
                raise FileNotFoundError(SCF_FOLDER_NOT_FOUND_ERROR.format(scf_folder_path))
            
            # Check if forward pass folders exist
            for folder_path in forward_pass_folder_path_list:
                if not os.path.exists(folder_path):
                    raise FileNotFoundError(FORWARD_PASS_FOLDER_NOT_FOUND_ERROR.format(folder_path))

            # Load quadrature nodes (radial grid points) and weights
            quadrature_nodes   = cls.load_quadrature_nodes_data(scf_folder_path)
            quadrature_weights = cls.load_quadrature_weights_data(scf_folder_path, len(quadrature_nodes))
            n_points = len(quadrature_nodes)

            # Load the features data for scf_functional
            rho = cls.load_rho_data(scf_folder_path, n_points)
            scf_feature_data = cls.load_features(
                folder_path                     = scf_folder_path,
                features_list                   = features_list,
                reference_array_length          = n_points,
                use_feature_round_off           = use_feature_round_off,
                smooth_reduced_laplacian        = smooth_reduced_laplacian,
                smooth_rho_for_derivatives      = smooth_rho_for_derivatives,
                smooth_rho_r_a                  = smooth_rho_r_a,
                smooth_rho_r_b                  = smooth_rho_r_b,
                smooth_rho_r_c                  = smooth_rho_r_c,
                smooth_rho_log_floor            = smooth_rho_log_floor,
            )

            # Load the xc potential data for scf_functional
            v_x = cls.load_v_x_data(scf_folder_path, n_points)
            v_c = cls.load_v_c_data(scf_folder_path, n_points)
            if include_energy_density:
                e_x = cls.load_e_x_data(scf_folder_path, n_points)
                e_c = cls.load_e_c_data(scf_folder_path, n_points)
            else:
                e_x = None
                e_c = None

            scf_xc_data : XCDataType = (v_x, v_c, e_x, e_c)

            if use_radius_cutoff:
                cutoff_idx = cls.compute_cutoff_index(
                    rho = rho, 
                    v_x = v_x, 
                    v_c = v_c, 
                    rho_threshold = radius_cutoff_rho_threshold, 
                    v_x_threshold = radius_cutoff_v_x_threshold, 
                    v_c_threshold = radius_cutoff_v_c_threshold
                )

                # Get cutoff radius
                cutoff_radius = quadrature_nodes[cutoff_idx - 1] if cutoff_idx > 0 else 0.0

                # Apply mask to features data, update n_points to filtered size
                n_filtered         = cutoff_idx
                scf_feature_data   = scf_feature_data[:cutoff_idx]
                if include_energy_density:
                    scf_xc_data = (v_x[:cutoff_idx], v_c[:cutoff_idx], e_x[:cutoff_idx], e_c[:cutoff_idx])
                else:
                    scf_xc_data = (v_x[:cutoff_idx], v_c[:cutoff_idx], None, None)
            else:
                # No cutoff filtering, use all data
                cutoff_radius = quadrature_nodes[-1]
                cutoff_idx    = n_points
                n_filtered    = n_points

            # Load the xc potential data for forward pass functionals
            forward_pass_xc_data_list : List[XCDataType] = []
            for forward_pass_folder_path in forward_pass_folder_path_list:
                v_x = cls.load_v_x_data(forward_pass_folder_path, n_points)
                v_c = cls.load_v_c_data(forward_pass_folder_path, n_points)
                if include_energy_density:
                    e_x = cls.load_e_x_data(forward_pass_folder_path, n_points)
                    e_c = cls.load_e_c_data(forward_pass_folder_path, n_points)
                    forward_pass_xc_data_list.append((v_x[:cutoff_idx], v_c[:cutoff_idx], e_x[:cutoff_idx], e_c[:cutoff_idx]))
                else:
                    forward_pass_xc_data_list.append((v_x[:cutoff_idx], v_c[:cutoff_idx], None, None))

            # Configuration ID (config_id) convention:
            #   0       = "main" folder (converged SCF result)
            #   1, 2,…  = "outer_iter_1", "outer_iter_2", … (intermediate iterations)
            #   -1      = unknown (non-standard folder label or outer_iter_XX parse failed)
            if folder_label == "main":
                config_id = 0
            elif folder_label.startswith("outer_iter_"):
                try:
                    config_id = int(folder_label.split("_")[-1])
                except ValueError:
                    config_id = -1  # Unknown intermediate
            else:
                config_id = -1  # Unknown label

            # Load derivative matrix if available (only for main configuration to save space)
            derivative_matrix_use_shared = False
            derivative_matrix = None
            if folder_label == "main":
                # Get configuration folder path (parent of scf_folder_path) and data_root
                config_folder = os.path.dirname(scf_folder_path)
                data_root = os.path.dirname(config_folder)  # Go up from configuration_XXX/atom_XXX to data_root
                
                # Use modular function to load derivative matrix
                derivative_matrix_use_shared, derivative_matrix = cls.load_derivative_matrix(
                    config_folder_path=config_folder,
                    data_root=data_root,
                )

            return SingleConfigurationData(
                atomic_number                   = atomic_number,
                n_electrons                     = n_electrons,
                configuration_id                = config_id,
                cutoff_radius                   = cutoff_radius,
                cutoff_idx                      = cutoff_idx,
                n_points                        = n_points,
                n_filtered                      = n_filtered,
                folder_label                    = folder_label,
                features_list                   = features_list,
                scf_xc_functional               = scf_xc_functional,
                forward_pass_xc_functional_list = forward_pass_xc_functional_list,
                scf_folder_path                 = scf_folder_path,
                forward_pass_folder_path_list   = forward_pass_folder_path_list,

                features_data                   = scf_feature_data,
                quadrature_nodes                = quadrature_nodes,
                quadrature_weights              = quadrature_weights,
                scf_xc_data                     = scf_xc_data,
                forward_pass_xc_data_list       = forward_pass_xc_data_list,
                derivative_matrix_use_shared    = derivative_matrix_use_shared,
                derivative_matrix               = derivative_matrix,
            )


        except Exception as e:
            # Import format_error_message locally to avoid circular import
            from .data_manager import format_error_message
            error_summary, error_traceback = format_error_message(e, f"Error loading data from {folder_label} folder")
            print(error_summary)
            print(f"Traceback:\n{error_traceback}")
            return None




    @classmethod
    def load_data(cls,
        # Required arguments
        data_root                   : str,
        scf_xc_functional           : str,
        forward_pass_xc_functionals : List[str],
        features_list               : List[str],
        configuration_index_list    : List[int],

        # Control arguments
        use_radius_cutoff           : bool = False,
        use_feature_round_off       : bool = True,
        include_energy_density      : bool = False,
        include_intermediate        : bool = False,
        print_debug_info            : bool = False,

        # Additional arguments, for parameter 'tuning'
        radius_cutoff_rho_threshold : Optional[float] = 1e-6,
        radius_cutoff_v_x_threshold : Optional[float] = 1e-8,
        radius_cutoff_v_c_threshold : Optional[float] = 1e-8,
        smooth_reduced_laplacian    : Optional[bool]  = False,
        smooth_rho_for_derivatives  : bool = False,
        smooth_rho_r_a              : float = 7.0,
        smooth_rho_r_b              : float = 10.0,
        smooth_rho_r_c              : float = 13.0,
        smooth_rho_log_floor        : float = 1e-30,
    ) -> Tuple[List[SingleConfigurationData], List[int]]:
        """
        Load training data for all configurations or specified ones.
        
        Parameters
        ----------
        data_root                   : Root directory of the dataset
        scf_xc_functional           : SCF XC functional (used for full SCF calculation to convergence)
        forward_pass_xc_functionals : Forward pass XC functional(s), if not None, will perform forward pass for each functional based on SCF results
        features_list               : List of features to load
        configuration_index_list    : List of configuration indices (1-based) to load. Atomic numbers are read from meta.json files.
        use_radius_cutoff           : Whether to apply cutoff filtering based on rho and vxc thresholds
        use_feature_round_off       : Whether to apply feature round-off (lower/upper bounds)
        include_energy_density      : Whether to also load energy density data
        include_intermediate        : Whether to also load data from intermediate iteration folders (outer_iter_XX)
        print_debug_info            : Whether to print debug information
        radius_cutoff_rho_threshold : Threshold for rho data when use_radius_cutoff is True
        radius_cutoff_v_x_threshold : Threshold for v_x data when use_radius_cutoff is True
        radius_cutoff_v_c_threshold : Threshold for v_c data when use_radius_cutoff is True
        smooth_reduced_laplacian    : Whether to smooth the reduced Laplacian at large radius using the weak-form 
        smooth_rho_for_derivatives  : If True, log-linear asymptotic rho smoothing then FEM-based grad/lap (see load_features);
            configurations whose smoothed rho is not monotone non-increasing on the radial grid are dropped with a warning.

        Returns
        -------
        configuration_data_list : List[SingleConfigurationData]
            List of SingleConfigurationData objects, 
            If the intermediate configurations are not included, the list will only contain one converged configuration data for each atom
        skipped_atoms : List[int]
            List of skipped atomic numbers
        """
        
        # Load data for each configuration
        skipped_atoms           : List[int]                     = []
        configuration_data_list : List[SingleConfigurationData] = []

        for config_idx in configuration_index_list:

            # Find configuration folder (with backward compatibility)
            config_folder = DataLoader._find_configuration_folder(data_root, config_idx)
            if config_folder is None:
                print(CONFIGURATION_FOLDER_NOT_FOUND_WARNING.format(config_idx))
                skipped_atoms.append(config_idx)
                continue
            
            scf_folder = os.path.join(config_folder, scf_xc_functional.lower())
            forward_pass_folder_path_list = [os.path.join(config_folder, fp_func.lower()) 
                for fp_func in forward_pass_xc_functionals]
            
            # Read atomic number, n_electrons, and converged status from meta.json
            atomic_number, n_electrons, converged = \
                DataLoader._load_metadata_from_meta(config_folder, config_idx)
            
            # Check convergence status - skip if not converged
            if converged is False:
                print(CONFIGURATION_NOT_CONVERGED_WARNING.format(config_idx, atomic_number))
                skipped_atoms.append(config_idx)
                continue
            
            # Check if SCF and forward pass folders exist
            if not os.path.exists(scf_folder):
                print(SCF_FOLDER_NOT_FOUND_WARNING.format(scf_folder, atomic_number))
                skipped_atoms.append(config_idx)
                continue
            
            for forward_pass_folder_path in forward_pass_folder_path_list:
                if not os.path.exists(forward_pass_folder_path):
                    print(FORWARD_PASS_FOLDER_NOT_FOUND_WARNING.format(forward_pass_folder_path, atomic_number))
                    skipped_atoms.append(config_idx)
                    continue

            # Load main folder data (forward_pass_folders can be empty list if forward_pass_xc_functional_list is empty)
            converged_configuration_data = cls.load_single_configuration_data(
                atomic_number                   = atomic_number,
                n_electrons                     = n_electrons,
                features_list                   = features_list,
                scf_xc_functional               = scf_xc_functional,
                forward_pass_xc_functional_list = forward_pass_xc_functionals,
                scf_folder_path                 = scf_folder,
                forward_pass_folder_path_list   = forward_pass_folder_path_list,
                use_radius_cutoff               = use_radius_cutoff,
                use_feature_round_off           = use_feature_round_off,
                include_energy_density          = include_energy_density,
                folder_label                    = "main",
                radius_cutoff_rho_threshold     = radius_cutoff_rho_threshold,
                radius_cutoff_v_x_threshold     = radius_cutoff_v_x_threshold,
                radius_cutoff_v_c_threshold     = radius_cutoff_v_c_threshold,
                smooth_reduced_laplacian        = smooth_reduced_laplacian,
                smooth_rho_for_derivatives      = smooth_rho_for_derivatives,
                smooth_rho_r_a                  = smooth_rho_r_a,
                smooth_rho_r_b                  = smooth_rho_r_b,
                smooth_rho_r_c                  = smooth_rho_r_c,
                smooth_rho_log_floor            = smooth_rho_log_floor,
            )
            
            if converged_configuration_data is None:
                print(CONVERGED_CONFIGURATION_DATA_NOT_FOUND_WARNING.format(atomic_number))
                skipped_atoms.append(config_idx)
                continue

            if smooth_rho_for_derivatives:
                rho_col = converged_configuration_data.features_data[
                    :, converged_configuration_data.features_list.index("rho")
                ]
                if not cls._rho_is_monotone_nonincreasing(rho_col):
                    msg = SMOOTH_RHO_NON_MONOTONE_DROPPED_WARNING.format(
                        config_idx, "main", atomic_number
                    )
                    warnings.warn(msg, UserWarning, stacklevel=2)
                    skipped_atoms.append(config_idx)
                    continue

            # Store converged configuration data
            configuration_data_list.append(converged_configuration_data)
            if print_debug_info:
                print(CONVERGED_CONFIGURATION_DATA_LOADED_DEBUG.format(
                    config_idx,
                    atomic_number, 
                    converged_configuration_data.n_electrons,
                    converged_configuration_data.n_points, 
                    converged_configuration_data.n_filtered, 
                    converged_configuration_data.cutoff_radius,
                ))
            
            # Load intermediate folders if requested
            if include_intermediate:
                # Find all outer_iter_XX folders in scf folder
                intermediate_folders = []
                if os.path.exists(scf_folder):
                    for item in os.listdir(scf_folder):
                        if item.startswith('outer_iter_') and os.path.isdir(os.path.join(scf_folder, item)):
                            intermediate_folders.append(item)
                
                # Sort intermediate folders by iteration number
                intermediate_folders.sort()
                
                if len(intermediate_folders) > 0:
                    if print_debug_info:
                        print(FOUND_INTERMEDIATE_CONFIGURATION_DATA_DEBUG.format(
                            len(intermediate_folders),
                            atomic_number,
                        ))
                    
                    for intermediate_folder in intermediate_folders:
                        intermediate_scf_path = os.path.join(scf_folder, intermediate_folder)
                        intermediate_forward_pass_folder_path_list = [os.path.join(fp_folder_path, intermediate_folder) 
                            for fp_folder_path in forward_pass_folder_path_list]
                        
                        # Load intermediate folder data (intermediate_forward_pass_paths can be empty list)
                        intermediate_configuration_data = cls.load_single_configuration_data(
                            atomic_number                   = atomic_number,
                            n_electrons                     = n_electrons,
                            features_list                   = features_list,
                            scf_xc_functional               = scf_xc_functional,
                            forward_pass_xc_functional_list = forward_pass_xc_functionals,
                            scf_folder_path                 = intermediate_scf_path,
                            forward_pass_folder_path_list   = intermediate_forward_pass_folder_path_list,
                            use_radius_cutoff               = use_radius_cutoff,
                            use_feature_round_off           = use_feature_round_off,
                            include_energy_density          = include_energy_density,
                            folder_label                    = intermediate_folder,
                            radius_cutoff_rho_threshold     = radius_cutoff_rho_threshold,
                            radius_cutoff_v_x_threshold     = radius_cutoff_v_x_threshold,
                            radius_cutoff_v_c_threshold     = radius_cutoff_v_c_threshold,
                            smooth_reduced_laplacian        = smooth_reduced_laplacian,
                            smooth_rho_for_derivatives      = smooth_rho_for_derivatives,
                            smooth_rho_r_a                  = smooth_rho_r_a,
                            smooth_rho_r_b                  = smooth_rho_r_b,
                            smooth_rho_r_c                  = smooth_rho_r_c,
                            smooth_rho_log_floor            = smooth_rho_log_floor,
                        )
                        
                        if intermediate_configuration_data is not None:
                            if smooth_rho_for_derivatives:
                                rho_col = intermediate_configuration_data.features_data[
                                    :,
                                    intermediate_configuration_data.features_list.index("rho"),
                                ]
                                if not cls._rho_is_monotone_nonincreasing(rho_col):
                                    msg = SMOOTH_RHO_NON_MONOTONE_DROPPED_WARNING.format(
                                        config_idx,
                                        intermediate_folder,
                                        atomic_number,
                                    )
                                    warnings.warn(msg, UserWarning, stacklevel=2)
                                    skipped_atoms.append(config_idx)
                                    continue
                            
                            # Convert data format for backward compatibility
                            configuration_data_list.append(intermediate_configuration_data)
                            if print_debug_info:
                                print(INTERMEDIATE_CONFIGURATION_DATA_LOADED_DEBUG.format(
                                    atomic_number,
                                    intermediate_configuration_data.n_points,
                                    intermediate_configuration_data.n_filtered,
                                    intermediate_configuration_data.n_electrons,
                                    intermediate_configuration_data.cutoff_radius,
                                ))
                        else:
                            print(INTERMEDIATE_CONFIGURATION_DATA_NOT_FOUND_WARNING.format(
                                atomic_number,
                                intermediate_folder,
                            ))
                            skipped_atoms.append(config_idx)
                            continue
        
        # Concatenate all data
        if len(configuration_data_list) == 0:
            raise ValueError(NO_DATA_LOADED_ERROR)
        
        # Print summary of skipped atoms if any
        if skipped_atoms:
            print(SKIPPED_ATOMS_WARNING.format(len(skipped_atoms), skipped_atoms))
        
        return configuration_data_list, skipped_atoms


