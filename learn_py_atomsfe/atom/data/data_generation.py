
__author__ = "Qihao Cheng"

"""Data generation utilities for atomic DFT calculations."""
import json
import numpy as np
import os
import sys
import traceback
from io import StringIO
from typing import Tuple, TextIO, TYPE_CHECKING
from pathlib import Path
from typing import Optional

if TYPE_CHECKING:
    from ..solver import AtomicDFTSolver

from ..scf.density import DensityCalculator
from ..mesh.operators import GridData

# Error messages
ATOMIC_NUMBER_NOT_INT_OR_FLOAT_ERROR = \
    "parameter 'atomic_number' must be an integer or float, get {} instead."
ATOMIC_NUMBER_NOT_INTEGER_WHEN_USED_AS_INDEX_ERROR = \
    "when 'configuration_index' is None, 'atomic_number' must be an integer to be used as configuration folder index, but got {} (non-integer float)."
N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR = \
    "parameter 'n_electrons' must be an integer or float, get {} instead."
N_ELECTRONS_IS_NONE_ERROR = \
    "parameter 'n_electrons' must not be None. It is a required parameter."
N_ELECTRONS_LIST_LENGTH_MISMATCH_ERROR = \
    "parameter 'n_electrons_list' length {} must match atomic_number_list length {}."
DOMAIN_SIZE_NOT_FLOAT_ERROR = \
    "parameter 'domain_size' must be a float, get {} instead."
FINITE_ELEMENTS_NUMBER_NOT_INTEGER_ERROR = \
    "parameter 'finite_elements_number' must be an integer, get {} instead."
POLYNOMIAL_ORDER_NOT_INTEGER_ERROR = \
    "parameter 'polynomial_order' must be an integer, get {} instead."
QUADRATURE_POINT_NUMBER_NOT_INTEGER_ERROR = \
    "parameter 'quadrature_point_number' must be an integer, get {} instead."
OEP_BASIS_NUMBER_NOT_INTEGER_ERROR = \
    "parameter 'oep_basis_number' must be an integer, get {} instead."
MAX_SCF_ITERATIONS_NOT_INTEGER_ERROR = \
    "parameter 'max_scf_iterations' must be an integer, get {} instead."
MAX_SCF_ITERATIONS_OUTER_NOT_INTEGER_ERROR = \
    "parameter 'max_scf_iterations_outer' must be an integer, get {} instead."
START_CONFIGURATION_INDEX_NOT_INTEGER_ERROR = \
    "parameter 'start_configuration_index' must be a positive integer, get {} instead."
START_CONFIGURATION_INDEX_NOT_POSITIVE_ERROR = \
    "parameter 'start_configuration_index' must be >= 1, get {} instead."
VERBOSE_NOT_BOOL_ERROR = \
    "parameter 'verbose' must be a boolean, get {} instead."
XC_FUNCTIONAL_NOT_STRING_ERROR = \
    "parameter 'xc_functional' must be a string, get {} instead."
SCF_XC_FUNCTIONAL_NOT_STRING_ERROR = \
    "parameter 'scf_xc_functional' must be a string, get {} instead."
FORWARD_PASS_XC_FUNCTIONAL_NOT_STRING_ERROR = \
    "parameter 'forward_pass_xc_functional' must be a string, get {} instead."
FORWARD_PASS_XC_FUNCTIONAL_NOT_STRING_OR_LIST_OF_STRINGS_ERROR = \
    "parameter 'forward_pass_xc_functionals' must be a string or list of strings, get type {} instead."
MESH_TYPE_NOT_STRING_ERROR = \
    "parameter 'mesh_type' must be a string, get {} instead."
MESH_CONCENTRATION_NOT_FLOAT_ERROR = \
    "parameter 'mesh_concentration' must be a float, get {} instead."
DIRECTORY_PATH_NOT_STRING_ERROR = \
    "parameter 'directory_path' must be a string, get {} instead."
SAVE_ENERGY_DENSITY_NOT_BOOL_ERROR = \
    "parameter 'save_energy_density' must be a boolean, get {} instead."
SAVE_INTERMEDIATE_NOT_BOOL_ERROR = \
    "parameter 'save_intermediate' must be a boolean, get {} instead."
PROCESS_INTERMEDIATE_NOT_BOOL_ERROR = \
    "parameter 'process_intermediate' must be a boolean, get {} instead."
SAVE_FULL_SPECTRUM_NOT_BOOL_ERROR = \
    "parameter 'save_full_spectrum' must be a boolean, get {} instead."
SAVE_FULL_ORBITALS_NOT_BOOL_ERROR = \
    "parameter 'save_full_orbitals' must be a boolean, get {} instead."

SAVE_INTERMEDIATE_INFO_IS_NONE_ERROR = \
    "intermediate_info should not be None when save_intermediate=True, get {} instead.\n \
    Remark: This error should never be raised, otherwise you should check the output of atomic_dft_solver.solve()."


# Warning messages
# Note: Use f-string or .format() with keyword arguments to avoid KeyError
OEP_IS_ENABLED_BUT_FULL_SPECTRUM_DATA_FILES_NOT_FOUND_WARNING = \
    "Warning: OEP is enabled but full spectrum data files not found in '{}'. Disabling OEP."

# Type aliases
# (v_x_local, v_c_local, e_x_local, e_c_local), where e_x_local and e_c_local are only available if compute_energy_density is True
XCDataType = Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]



class TeeOutput:
    """Custom output class that writes to both terminal and buffer simultaneously.
    
    When verbose=False, only writes to buffer (for saving to out.txt), not to terminal.
    When verbose=True, writes to both terminal and buffer.
    """
    def __init__(self, terminal, buffer, verbose=True):
        self.terminal = terminal
        self.buffer = buffer
        self.verbose = verbose
    
    def write(self, message):
        # Always write to buffer (for out.txt file)
        self.buffer.write(message)
        # Only write to terminal if verbose is True
        if self.verbose:
            self.terminal.write(message)
    
    def flush(self):
        self.buffer.flush()
        if self.verbose:
            self.terminal.flush()


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


# Maps XC functionals to their default use_oep values (True or False)
# Principle: use OEP if the functional supports it
XC_FUNCTIONAL_OEP_DEFAULT = {
    "None"    : False,  # No XC functional
    "LDA_PZ"  : False,  # LDA Perdew-Zunger
    "LDA_PW"  : False,  # LDA Perdew-Wang
    "GGA_PBE" : False,  # GGA Perdew-Burke-Ernzerhof
    "SCAN"    : False,  # SCAN functional, meta-GGA
    "RSCAN"   : False,  # RSCAN functional, meta-GGA
    "R2SCAN"  : False,  # R2SCAN functional, meta-GGA
    "HF"      : False,  # Hartree-Fock
    "PBE0"    : True,   # PBE0 Perdew-Burke-Ernzerhof, hybrid functional (supports OEP)
    "EXX"     : True,   # Exact Exchange, using OEP method
    "RPA"     : True,   # Random Phase Approximation, with exact exchange
}


class DataGenerator:
    @classmethod
    def generate_single_atom_data(
        cls,
        # Required arguments
        atomic_number            : int | float,
        n_electrons              : Optional[int | float],
        use_oep                  : bool,
        directory_path           : str,
        configuration_index      : int,

        # Arguments controlling the contents of the dataset
        save_energy_density      : bool,
        save_intermediate        : bool,
        save_full_spectrum       : bool,
        save_full_orbitals       : bool,
        save_derivative_matrix   : bool,

        # Arguments controlling the generation process
        # Grid, basis, and mesh parameters
        domain_size              : float,
        finite_elements_number   : int,
        polynomial_order         : int,
        quadrature_point_number  : int,
        oep_basis_number         : int,
        mesh_type                : str,
        mesh_concentration       : float,
        mesh_spacing             : float,

        # SCF convergence parameters
        xc_functional            : str,
        scf_tolerance            : float,
        max_scf_iterations       : int,
        max_scf_iterations_outer : int,
        use_pulay_mixing         : bool,
        use_preconditioner       : bool,
        pulay_mixing_parameter   : float,
        pulay_mixing_history     : int,
        pulay_mixing_frequency   : int,
        linear_mixing_alpha1     : float,
        linear_mixing_alpha2     : float,

        # Advanced functional parameters
        hybrid_mixing_parameter           : float,
        frequency_quadrature_point_number : int,
        angular_momentum_cutoff           : int,
        double_hybrid_flag                : bool,
        oep_mixing_parameter              : float,
        enable_parallelization            : bool,

        # Debugging and verbose parameters
        verbose                  : bool,

        # Derivative matrix handling
        shared_derivative_matrix_path : Optional[str] = None,
    ) -> None:
        """
        Generate data for single atomic number.
        
        Required arguments
        ------------------
        `atomic_number` : int | float
            Atomic number.
        `n_electrons` : int | float | None
            Number of electrons. If None, defaults to atomic_number.
        `use_oep` : bool
            Whether to use OEP for the SCF calculation.
        `directory_path` : str
            Directory path for saving data.
        `configuration_index` : int
            Index of the configuration in the list.

        Arguments controlling the contents of the dataset
        -------------------------------------------------
        `save_energy_density` : bool
            Whether to save energy density.
        `save_intermediate` : bool
            Whether to save intermediate information.
        `save_full_spectrum` : bool
            Whether to save full_eigen_energies and full_l_terms.
        `save_full_orbitals` : bool
            Whether to save full_orbitals (default: False, to save storage). Only has effect when full spectrum is computed.
        `save_derivative_matrix` : bool
            Whether to save derivative matrix. Most systems have the same derivative matrix when using
            the same grid/basis/mesh parameters, so a shared derivative matrix is saved at the dataset root.
            If an atom's derivative matrix differs from the shared one, it is saved locally and recorded in meta.json.

        Arguments controlling the generation process
        ---------------------------------------------
        Grid, basis, and mesh parameters
        `domain_size` : float
            Radial computational domain size in atomic units.
        `finite_elements_number` : int
            Number of finite elements.
        `polynomial_order` : int
            Polynomial order of basis functions.
        `quadrature_point_number` : int
            Number of quadrature points.
        `oep_basis_number` : int
            OEP basis number.
        `mesh_type` : str
            Mesh distribution type.
        `mesh_concentration` : float
            Mesh concentration parameter.
        `mesh_spacing` : float
            Output uniform mesh spacing.

        SCF convergence parameters
        --------------------------
        `xc_functional` : str
            XC functional.
        `scf_tolerance` : float
            SCF convergence tolerance.
        `max_scf_iterations` : int | None
            Maximum number of inner SCF iterations. If None, uses default (500).
        `max_scf_iterations_outer` : int | None
            Maximum number of outer SCF iterations (for functionals requiring outer loop like HF, EXX, RPA, PBE0).
            If None, uses default (50 when needed, otherwise not used).
        `use_pulay_mixing` : bool
            True for Pulay mixing, False for linear mixing.
        `use_preconditioner` : bool
            Whether to use preconditioner.
        `pulay_mixing_parameter` : float
            Pulay mixing parameter.
        `pulay_mixing_history` : int
            Pulay mixing history.
        `pulay_mixing_frequency` : int
            Pulay mixing frequency.
        `linear_mixing_alpha1` : float
            Linear mixing parameter (alpha_1).
        `linear_mixing_alpha2` : float
            Linear mixing parameter (alpha_2).

        Debugging and verbose parameters
        --------------------------------
        `verbose` : bool
            Whether to print debug information.
        """

        # Type checks only, other checks are done in the AtomicDFTSolver class
        if not isinstance(atomic_number, (int, float)):
            raise TypeError(ATOMIC_NUMBER_NOT_INT_OR_FLOAT_ERROR.format(atomic_number))
        if n_electrons is None:
            n_electrons = atomic_number
        if not isinstance(n_electrons, (int, float)):
            raise TypeError(N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(n_electrons))
        n_electrons = float(n_electrons)
        if not isinstance(domain_size, float):
            raise TypeError(DOMAIN_SIZE_NOT_FLOAT_ERROR.format(domain_size))
        if not isinstance(finite_elements_number, int):
            raise TypeError(FINITE_ELEMENTS_NUMBER_NOT_INTEGER_ERROR.format(finite_elements_number))
        if not isinstance(polynomial_order, int):
            raise TypeError(POLYNOMIAL_ORDER_NOT_INTEGER_ERROR.format(polynomial_order))
        if not isinstance(quadrature_point_number, int):
            raise TypeError(QUADRATURE_POINT_NUMBER_NOT_INTEGER_ERROR.format(quadrature_point_number))
        if not isinstance(oep_basis_number, int):
            raise TypeError(OEP_BASIS_NUMBER_NOT_INTEGER_ERROR.format(oep_basis_number))
        if max_scf_iterations is not None and not isinstance(max_scf_iterations, int):
            raise TypeError(MAX_SCF_ITERATIONS_NOT_INTEGER_ERROR.format(type(max_scf_iterations)))
        if max_scf_iterations_outer is not None and not isinstance(max_scf_iterations_outer, int):
            raise TypeError(MAX_SCF_ITERATIONS_OUTER_NOT_INTEGER_ERROR.format(type(max_scf_iterations_outer)))
        if not isinstance(verbose, bool):
            raise TypeError(VERBOSE_NOT_BOOL_ERROR.format(verbose))
        if not isinstance(xc_functional, str):
            raise TypeError(XC_FUNCTIONAL_NOT_STRING_ERROR.format(xc_functional))
        if not isinstance(directory_path, str):
            raise TypeError(DIRECTORY_PATH_NOT_STRING_ERROR.format(directory_path))
        if not isinstance(mesh_type, str):
            raise TypeError(MESH_TYPE_NOT_STRING_ERROR.format(mesh_type))
        if not isinstance(mesh_concentration, float):
            raise TypeError(MESH_CONCENTRATION_NOT_FLOAT_ERROR.format(mesh_concentration))
        if not isinstance(save_energy_density, bool):
            raise TypeError(SAVE_ENERGY_DENSITY_NOT_BOOL_ERROR.format(save_energy_density))
        if not isinstance(save_intermediate, bool):
            raise TypeError(SAVE_INTERMEDIATE_NOT_BOOL_ERROR.format(save_intermediate))
        if not isinstance(save_full_spectrum, bool):
            raise TypeError(SAVE_FULL_SPECTRUM_NOT_BOOL_ERROR.format(save_full_spectrum))
        if not isinstance(save_full_orbitals, bool):
            raise TypeError(SAVE_FULL_ORBITALS_NOT_BOOL_ERROR.format(save_full_orbitals))

        # Create folder for this configuration (index-based, not atomic number)
        # New format: configuration_XXX (replaces old atom_XXX format)
        config_folder = os.path.join(directory_path, f"configuration_{configuration_index:03d}")
        folder_name = os.path.join(config_folder, xc_functional.lower())
        os.makedirs(folder_name, exist_ok=True)

        # Save minimal metadata at configuration folder level
        os.makedirs(config_folder, exist_ok=True)
        meta_path = os.path.join(config_folder, "meta.json")
        with open(meta_path, "w") as meta_file:
            json.dump(
                {"atomic_number": float(atomic_number), "n_electrons": n_electrons},
                meta_file,
                indent=2,
                sort_keys=True
            )

        # Setup output redirection to both terminal and buffer
        # When verbose=False, only capture to buffer (for out.txt), not to terminal
        output_buffer = StringIO()
        original_stdout = sys.stdout
        tee_output = TeeOutput(original_stdout, output_buffer, verbose=verbose)
        sys.stdout = tee_output

        try:
            # Run atomic DFT calculation
            # Always use verbose=True for AtomicDFTSolver to ensure detailed output is captured in buffer
            # TeeOutput will control whether output goes to terminal based on verbose parameter
            from ..solver import AtomicDFTSolver
                        
            atomic_dft_solver = AtomicDFTSolver(
                # Required arguments
                atomic_number             = atomic_number,
                n_electrons               = n_electrons,
                all_electron_flag         = True, # Always do all-electron calculations
                xc_functional             = xc_functional,
                use_oep                   = use_oep,

                # Grid, basis, and mesh parameters
                domain_size               = domain_size,
                finite_element_number     = finite_elements_number,
                polynomial_order          = polynomial_order,
                quadrature_point_number   = quadrature_point_number,
                oep_basis_number          = oep_basis_number,
                mesh_type                 = mesh_type,
                mesh_concentration        = mesh_concentration,
                mesh_spacing              = mesh_spacing,

                # SCF convergence parameters
                scf_tolerance             = scf_tolerance,
                max_scf_iterations        = max_scf_iterations,
                max_scf_iterations_outer  = max_scf_iterations_outer,
                use_pulay_mixing          = use_pulay_mixing,
                use_preconditioner        = use_preconditioner,
                pulay_mixing_parameter    = pulay_mixing_parameter,
                pulay_mixing_history      = pulay_mixing_history,
                pulay_mixing_frequency    = pulay_mixing_frequency,
                linear_mixing_alpha1      = linear_mixing_alpha1,
                linear_mixing_alpha2      = linear_mixing_alpha2,

                # Advanced functional parameters
                hybrid_mixing_parameter           = hybrid_mixing_parameter,
                frequency_quadrature_point_number = frequency_quadrature_point_number,
                angular_momentum_cutoff           = angular_momentum_cutoff,
                double_hybrid_flag                = double_hybrid_flag,
                oep_mixing_parameter              = oep_mixing_parameter,
                enable_parallelization            = enable_parallelization,

                # Debugging and verbose parameters
                verbose                   = True,  # Always True to capture detailed output in buffer 
            )
            final_result = atomic_dft_solver.solve(
                save_intermediate   = save_intermediate, 
                save_energy_density = save_energy_density,
                save_full_spectrum  = save_full_spectrum or save_full_orbitals,
                evaluate_basis_on_uniform_grid = True,
            )
            
            # Update meta.json with convergence information
            converged = final_result.get('converged', False)
            with open(meta_path, "r") as meta_file:
                meta_data = json.load(meta_file)
            meta_data['converged'] = bool(converged)
            with open(meta_path, "w") as meta_file:
                json.dump(meta_data, meta_file, indent=2, sort_keys=True)
            
            # Restore stdout
            sys.stdout = original_stdout
            
            # Save captured output to out.txt
            with open(os.path.join(folder_name, "out.txt"), "w") as f:
                f.write(output_buffer.getvalue())
            
            # Helper function to save data to a folder, use existing solver instance to avoid creating new one
            def save_data_to_folder(target_folder, rho_data, orbitals_data, v_x_data, v_c_data, 
                                    quadrature_nodes_data, quadrature_nodes_uniform_data, quadrature_weights_data,
                                    v_x_uniform_data, v_c_uniform_data):
                """Save all data files to the specified folder."""
                grad_rho = atomic_dft_solver.scf_driver.density_calculator.compute_density_gradient(rho_data)
                lap_rho = atomic_dft_solver.scf_driver.density_calculator.compute_density_laplacian(rho_data)
                hartree = atomic_dft_solver.scf_driver.poisson_solver.solve_hartree(rho_data)
                
                np.savetxt(os.path.join(target_folder, "quadrature_nodes.txt"), quadrature_nodes_data)
                np.savetxt(os.path.join(target_folder, "quadrature_weights.txt"), quadrature_weights_data)
                np.savetxt(os.path.join(target_folder, "rho.txt"), rho_data)
                np.savetxt(os.path.join(target_folder, "grad_rho.txt"), grad_rho)
                np.savetxt(os.path.join(target_folder, "lap_rho.txt"), lap_rho)
                np.savetxt(os.path.join(target_folder, "v_x.txt"), v_x_data)
                np.savetxt(os.path.join(target_folder, "v_c.txt"), v_c_data)
                np.savetxt(os.path.join(target_folder, "hartree.txt"), hartree)
                np.savetxt(os.path.join(target_folder, "orbitals.txt"), orbitals_data)
                
                np.savetxt(os.path.join(target_folder, "quadrature_nodes_uniform.txt"), quadrature_nodes_uniform_data)
                np.savetxt(os.path.join(target_folder, "v_x_uniform.txt"), v_x_uniform_data)
                np.savetxt(os.path.join(target_folder, "v_c_uniform.txt"), v_c_uniform_data)
                
            # Extract data from final results (converged state)
            quadrature_nodes = atomic_dft_solver.grid_data_standard.quadrature_nodes
            quadrature_weights = atomic_dft_solver.grid_data_standard.quadrature_weights
            rho = final_result['rho']
            v_x = final_result['v_x_local']
            v_c = final_result['v_c_local']
            orbitals = final_result['orbitals']
            
            # Save the v_x_pbe0 and v_c_pbe0 on the uniform grid
            quadrature_nodes_uniform = final_result['uniform_grid']
            v_x_uniform = final_result['v_x_local_on_uniform_grid']
            v_c_uniform = final_result['v_c_local_on_uniform_grid']
                        
            # Save the data to the folder
            save_data_to_folder(
                folder_name, rho, orbitals, v_x, v_c,
                quadrature_nodes, quadrature_nodes_uniform, quadrature_weights,
                v_x_uniform, v_c_uniform
            )

            # Save full_eigen_energies, full_orbitals, full_l_terms if available and requested (needed for OEP)
            # These arrays are very large, so only save if explicitly requested
            if save_full_spectrum:
                if final_result.get('full_eigen_energies') is not None:
                    np.savetxt(os.path.join(folder_name, "full_eigen_energies.txt"), final_result['full_eigen_energies'])
                if final_result.get('full_l_terms') is not None:
                    np.savetxt(os.path.join(folder_name, "full_l_terms.txt"), final_result['full_l_terms'])

            if save_full_orbitals:
                if final_result.get('full_orbitals') is not None:
                    np.savetxt(os.path.join(folder_name, "full_orbitals.txt"), final_result['full_orbitals'])

            if save_energy_density:
                # Save local XC energy density
                np.savetxt(os.path.join(folder_name, "e_x.txt"), final_result['e_x_local'])
                np.savetxt(os.path.join(folder_name, "e_c.txt"), final_result['e_c_local'])
                # Save local XC energy density on uniform grid
                np.savetxt(os.path.join(folder_name, "e_x_uniform.txt"), final_result['e_x_local_on_uniform_grid'])
                np.savetxt(os.path.join(folder_name, "e_c_uniform.txt"), final_result['e_c_local_on_uniform_grid'])

            # Handle derivative matrix saving if requested
            if save_derivative_matrix:
                cls._save_derivative_matrix(
                    derivative_matrix             = atomic_dft_solver.ops_builder_standard.get_derivative_matrix(),
                    shared_derivative_matrix_path = shared_derivative_matrix_path,
                    folder_name                   = folder_name,
                    directory_path                = directory_path,
                    meta_path                     = meta_path,
                )


            # save intermediate information if requested
            if save_intermediate:
                assert final_result['intermediate_info'] is not None, \
                    SAVE_INTERMEDIATE_INFO_IS_NONE_ERROR.format(final_result['intermediate_info'])
                intermediate_info = final_result['intermediate_info']

                # Save each outer iteration to a separate subfolder
                for outer_iter in intermediate_info.outer_iterations:
                    outer_iter_folder = os.path.join(folder_name, f"outer_iter_{outer_iter.outer_iteration:02d}")
                    os.makedirs(outer_iter_folder, exist_ok=True)

                    # Extract data from outer iteration
                    outer_rho          = outer_iter.density_data.rho
                    outer_orbitals     = outer_iter.orbitals
                    outer_density_data = outer_iter.density_data

                    # Compute XC potentials for this outer iteration using energy_calculator
                    # This method handles both local XC and OEP potentials correctly
                    outer_v_x, outer_v_c = atomic_dft_solver.energy_calculator.compute_local_xc_potential(
                        density_data           = outer_density_data,
                        full_eigen_energies    = outer_iter.full_eigen_energies,
                        full_orbitals          = outer_iter.full_orbitals,
                        full_l_terms           = outer_iter.full_l_terms,
                        enable_parallelization = atomic_dft_solver.enable_parallelization,
                    )

                    # Evaluate on uniform grid
                    outer_v_x_uniform = atomic_dft_solver.ops_builder_standard.evaluate_single_field_on_grid(
                        given_grid   = quadrature_nodes_uniform,
                        field_values = outer_v_x,
                    )
                    outer_v_c_uniform = atomic_dft_solver.ops_builder_standard.evaluate_single_field_on_grid(
                        given_grid   = quadrature_nodes_uniform,
                        field_values = outer_v_c,
                    )

                    if save_energy_density:
                        # Compute local XC energy density
                        outer_e_x, outer_e_c = atomic_dft_solver.energy_calculator.compute_local_xc_energy_density(
                            density_data           = outer_density_data,
                            full_eigen_energies    = outer_iter.full_eigen_energies,
                            full_orbitals          = outer_iter.full_orbitals,
                            full_l_terms           = outer_iter.full_l_terms,
                            enable_parallelization = atomic_dft_solver.enable_parallelization,
                        )

                        # Evaluate on uniform grid
                        outer_e_x_uniform = atomic_dft_solver.ops_builder_standard.evaluate_single_field_on_grid(
                            given_grid   = quadrature_nodes_uniform,
                            field_values = outer_e_x,
                        )
                        outer_e_c_uniform = atomic_dft_solver.ops_builder_standard.evaluate_single_field_on_grid(
                            given_grid   = quadrature_nodes_uniform,
                            field_values = outer_e_c,
                        )
                    else:
                        outer_e_x, outer_e_c = None, None
                        outer_e_x_uniform, outer_e_c_uniform = None, None

                    # Save the data to the folder
                    save_data_to_folder(
                        outer_iter_folder, outer_rho, outer_orbitals, outer_v_x, outer_v_c,
                        quadrature_nodes, quadrature_nodes_uniform, quadrature_weights,
                        outer_v_x_uniform, outer_v_c_uniform
                    )

                    # Save full_eigen_energies / full_l_terms (save_full_spectrum) and full_orbitals (save_full_orbitals)
                    if save_full_spectrum:
                        if outer_iter.full_eigen_energies is not None:
                            np.savetxt(os.path.join(outer_iter_folder, "full_eigen_energies.txt"), outer_iter.full_eigen_energies)
                        if outer_iter.full_l_terms is not None:
                            np.savetxt(os.path.join(outer_iter_folder, "full_l_terms.txt"), outer_iter.full_l_terms)
                    if save_full_orbitals and outer_iter.full_orbitals is not None:
                        np.savetxt(os.path.join(outer_iter_folder, "full_orbitals.txt"), outer_iter.full_orbitals)

                    if save_energy_density:
                        # Save local XC energy density
                        np.savetxt(os.path.join(outer_iter_folder, "e_x.txt"), outer_e_x)
                        np.savetxt(os.path.join(outer_iter_folder, "e_c.txt"), outer_e_c)
                        # Save local XC energy density on uniform grid
                        np.savetxt(os.path.join(outer_iter_folder, "e_x_uniform.txt"), outer_e_x_uniform)
                        np.savetxt(os.path.join(outer_iter_folder, "e_c_uniform.txt"), outer_e_c_uniform)

                    # Save additional metadata for this outer iteration
                    with open(os.path.join(outer_iter_folder, "metadata.txt"), "w") as f:
                        f.write(f"Outer iteration: {outer_iter.outer_iteration}\n")
                        f.write(f"Outer residual: {outer_iter.outer_rho_residual:.10e}\n")
                        f.write(f"Converged: {outer_iter.converged}\n")
                        f.write(f"Inner iterations: {outer_iter.iterations}\n")
                        f.write(f"Number of inner iterations recorded: {len(outer_iter.inner_iterations)}\n")
            print(f"Successfully generated data for atom {atomic_number} in folder '{folder_name}'")

        except Exception as e:
            # Save output log even on error
            os.makedirs(folder_name, exist_ok=True)
            current_output = output_buffer.getvalue()
            with open(os.path.join(folder_name, "out.txt"), "w") as f:
                f.write(current_output)
                f.write(f"\n\nError occurred: {str(e)}\n")
            raise


    @staticmethod
    def _save_derivative_matrix(
        derivative_matrix             : np.ndarray,
        shared_derivative_matrix_path : Optional[str],
        folder_name                   : str,
        directory_path                : str,
        meta_path                     : str,
    ) -> None:
        """
        Save derivative matrix with shared/unique logic.
        
        Most systems have the same derivative matrix when using the same grid/basis/mesh parameters,
        so we save a shared derivative matrix at the dataset root and only save unique ones per atom if different.
        
        Args:
            derivative_matrix             : Derivative matrix from ops_builder_standard (3D array)
            shared_derivative_matrix_path : Path to shared derivative matrix at dataset root (None if not saving)
            folder_name                   : Folder path for this configuration's functional (e.g., configuration_001/pbe0)
            directory_path                : Root directory of the dataset
            meta_path                     : Path to meta.json file for this configuration
        
        The function will:
        1. Compare current derivative matrix with shared one (if exists)
        2. If they match: mark as use_shared=True in meta.json, no local file saved
        3. If they differ: save unique derivative matrix locally and mark as use_shared=False in meta.json
        4. If shared matrix doesn't exist: save it as the shared one and mark current atom as use_shared=True
        """
        # Load or create shared derivative matrix
        use_shared = False
        if shared_derivative_matrix_path is not None:
            if os.path.exists(shared_derivative_matrix_path):
                # Load existing shared derivative matrix and compare with current one
                shared_derivative_matrix = np.load(shared_derivative_matrix_path)
                # Compare with current derivative matrix (using tight tolerance for numerical comparison)
                if np.allclose(derivative_matrix, shared_derivative_matrix, rtol=1e-10, atol=1e-12):
                    use_shared = True
            else:
                # First time: save current derivative matrix as the shared one at dataset root
                os.makedirs(os.path.dirname(shared_derivative_matrix_path), exist_ok=True)
                np.save(shared_derivative_matrix_path, derivative_matrix)
                use_shared = True
        
        # Update meta.json with derivative matrix information
        with open(meta_path, "r") as meta_file:
            meta_data = json.load(meta_file)
        
        if use_shared:
            # This atom uses the shared derivative matrix (no need to save a copy)
            meta_data["derivative_matrix"] = {
                "use_shared": True,
                "path": "derivative_matrix.npy"  # Relative to data_root
            }
        else:
            # This atom has a unique derivative matrix, save it locally
            unique_derivative_matrix_path = os.path.join(folder_name, "derivative_matrix.npy")
            np.save(unique_derivative_matrix_path, derivative_matrix)
            # Get relative path from data_root (directory_path) for cross-platform compatibility
            rel_path = os.path.relpath(unique_derivative_matrix_path, directory_path)
            meta_data["derivative_matrix"] = {
                "use_shared": False,
                "path": rel_path.replace("\\", "/")  # Use forward slashes for cross-platform compatibility
            }
        
        with open(meta_path, "w") as meta_file:
            json.dump(meta_data, meta_file, indent=2, sort_keys=True)


    @staticmethod
    def calculate_laplacian_in_weak_form(
        rho                : np.ndarray,
        solver             : "AtomicDFTSolver",
        sparser_grid_ratio : float = 0.5,
    ) -> np.ndarray:
        """
        Calculate Laplacian in weak form on a sparser grid.
            Typically used for calculating the reduced Laplacian,
            since smaller rho values are more sensitive to numerical errors.

        Parameters
        ----------
        rho : np.ndarray
            Electron density
        solver : AtomicDFTSolver
            AtomicDFT solver
        sparser_grid_ratio : float
            Ratio of the sparser grid to the original grid

        Returns
        -------
        np.ndarray
            Laplacian in weak form
        """
        original_grid_data          : GridData          = solver.grid_data_standard
        original_density_calculator : DensityCalculator = solver.density_calculator

        sparser_grid_data = GridData.from_basic(
            domain_size             = solver.domain_size,
            finite_element_number   = original_grid_data.finite_element_number,
            polynomial_order        = int(solver.polynomial_order * sparser_grid_ratio),
            quadrature_point_number = solver.quadrature_point_number,
            mesh_type               = solver.mesh_type,
            mesh_concentration      = solver.mesh_concentration,
        )

        sparser_density_calculator = DensityCalculator(
            grid_data         = sparser_grid_data,
            occupation_info   = original_density_calculator.occupation_info,
            derivative_matrix = original_density_calculator.derivative_matrix,
            verbose           = False,
        )

        laplacian_in_weak_form = sparser_density_calculator.compute_density_laplacian_in_weak_form(rho)

        return laplacian_in_weak_form


    @staticmethod
    def _forward_pass_single_folder(
        # Required arguments
        atomic_number           : int | float,
        n_electrons             : int | float,
        xc_functional           : str,
        read_folder_path        : str,
        write_folder_path       : str,

        # Arguments controlling the contents of the dataset
        compute_energy_density  : bool = False,
        save_full_spectrum      : bool = False,
        save_full_orbitals      : bool = False,

        # Arguments controlling the generation process
        # Grid, basis, and mesh parameters, but no SCF convergence parameters, because we are only doing forward pass
        domain_size             : float = 20.0,
        finite_elements_number  : int = 35,
        polynomial_order        : int = 20,
        quadrature_point_number : int = 43,
        oep_basis_number        : int = 5,
        mesh_type               : str = "polynomial",
        mesh_concentration      : float = 2.0,
        mesh_spacing            : float = 0.1,

        # Debugging and verbose parameters
        verbose                 : bool = True,

        # Optional parameters control
        output_buffer           : Optional[StringIO] = None,
        original_stdout         : Optional[TextIO] = None,
        atomic_dft_solver       : Optional["AtomicDFTSolver"] = None,
        save_output_log         : bool = True,
    ) -> Tuple[XCDataType, "AtomicDFTSolver"]:

        """
        Perform forward pass for a single folder.
        
        Parameters
        ----------
            atomic_number           : Atomic number
            n_electrons             : Number of electrons
            xc_functional           : XC functional
            read_folder_path        : Path to folder containing orbitals.txt
            write_folder_path       : Path to folder where results will be saved
            compute_energy_density  : Whether to compute energy density (default: False)
            save_full_spectrum      : If True, save full_eigen_energies and full_l_terms (default: False)
            save_full_orbitals     : If True, save full_orbitals (default: False, to save storage)
            domain_size             : Domain size (default: 20.0)
            finite_elements_number  : Number of finite elements (default: 35)
            polynomial_order        : Polynomial order (default: 20)
            quadrature_point_number : Number of quadrature points (default: 43)
            oep_basis_number        : OEP basis number (default: 5)
            mesh_type               : Mesh type (default: "polynomial")
            mesh_concentration      : Mesh concentration (default: 2.0)
            mesh_spacing            : Output uniform mesh spacing (default: 0.1)
            verbose                 : Whether to print information during execution (default: True)
            output_buffer           : Optional StringIO buffer for output capture (default: None)
            original_stdout         : Optional original stdout (for nested calls) (default: None)
            atomic_dft_solver       : Optional pre-created solver instance to reuse (default: None)
            save_output_log         : If True, save out.txt file in the write_folder_path (default: True)
        
        Returns
        -------
            Tuple[XCDataType, AtomicDFTSolver]: (v_x_local, v_c_local, e_x_local, e_c_local), atomic_dft_solver
        """

        # Type checks for required parameters
        if not isinstance(atomic_number, (int, float)):
            raise TypeError(ATOMIC_NUMBER_NOT_INT_OR_FLOAT_ERROR.format(atomic_number))
        if n_electrons is None:
            raise ValueError(N_ELECTRONS_IS_NONE_ERROR)
        if not isinstance(n_electrons, (int, float)):
            raise TypeError(N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(n_electrons))
        
        # This is a internal static method, so we don't need to check the type of other arguments
        # Setup output redirection if not already set
        # When verbose=False, only capture to buffer (for out.txt), not to terminal
        if output_buffer is None:
            output_buffer = StringIO()
            original_stdout = sys.stdout
            tee_output = TeeOutput(original_stdout, output_buffer, verbose=verbose)
            sys.stdout = tee_output
            should_restore_stdout = True
        else:
            should_restore_stdout = False
        
        try:
            
            # Load orbitals
            orbitals_file = os.path.join(read_folder_path, "orbitals.txt")
            if not os.path.exists(orbitals_file):
                raise FileNotFoundError(f"Orbitals file not found: {orbitals_file}")
            orbitals = np.loadtxt(orbitals_file)
            if orbitals.ndim == 1:
                orbitals = orbitals.reshape(-1, 1)
            
            # Get default use_oep value from dictionary based on xc_functional
            # Forward pass always uses default OEP value
            use_oep = XC_FUNCTIONAL_OEP_DEFAULT.get(xc_functional, False)
            
            # Load full_eigen_energies, full_orbitals, full_l_terms if they exist (only needed for RPA)
            # If files don't exist and use_oep=True, we need to disable OEP to avoid errors
            full_eigen_energies = None
            full_orbitals = None
            full_l_terms = None
            actual_use_oep = use_oep
            
            # Try to load full spectrum data if available (optional, only for RPA)
            full_eigen_energies_file = os.path.join(read_folder_path, "full_eigen_energies.txt")
            full_orbitals_file       = os.path.join(read_folder_path, "full_orbitals.txt")
            full_l_terms_file        = os.path.join(read_folder_path, "full_l_terms.txt")
            
            # Try to load full spectrum data if available (optional, only for RPA)
            if os.path.exists(full_eigen_energies_file):
                full_eigen_energies = np.loadtxt(full_eigen_energies_file)
            if os.path.exists(full_orbitals_file):
                full_orbitals = np.loadtxt(full_orbitals_file)
                if full_orbitals.ndim == 1:
                    full_orbitals = full_orbitals.reshape(-1, 1)
            if os.path.exists(full_l_terms_file):
                full_l_terms = np.loadtxt(full_l_terms_file)
                if full_l_terms.ndim == 1:
                    full_l_terms = full_l_terms.reshape(-1, 1)
                    
            # If OEP is enabled but files don't exist, disable OEP to avoid errors
            if actual_use_oep and (full_eigen_energies is None or full_orbitals is None or full_l_terms is None):
                actual_use_oep = False
                # Use position-based format to avoid KeyError issues
                print(OEP_IS_ENABLED_BUT_FULL_SPECTRUM_DATA_FILES_NOT_FOUND_WARNING.format(read_folder_path))
                
            # Create new solver instance if not provided
            if atomic_dft_solver is None or (atomic_dft_solver.use_oep != actual_use_oep):
                from ..solver import AtomicDFTSolver
                atomic_dft_solver = AtomicDFTSolver(
                    atomic_number             = atomic_number,
                    n_electrons               = n_electrons,
                    all_electron_flag         = True,
                    xc_functional             = xc_functional,
                    use_oep                   = actual_use_oep,

                    domain_size               = domain_size,
                    finite_element_number     = finite_elements_number,
                    polynomial_order          = polynomial_order,
                    quadrature_point_number   = quadrature_point_number,
                    oep_basis_number          = oep_basis_number,
                    mesh_type                 = mesh_type,
                    mesh_concentration        = mesh_concentration,
                    mesh_spacing              = mesh_spacing,

                    verbose                   = True,  # Always True to capture detailed output in buffer
                )

            # Perform forward pass
            final_result = atomic_dft_solver.forward(
                orbitals               = orbitals,
                full_eigen_energies    = full_eigen_energies,
                full_orbitals          = full_orbitals,
                full_l_terms           = full_l_terms,
                compute_energy_density = compute_energy_density,
            )
            v_x_local = final_result['v_x_local']
            v_c_local = final_result['v_c_local']
            e_x_local = final_result['e_x_local']
            e_c_local = final_result['e_c_local']

            # Get original grid nodes
            quadrature_nodes   = atomic_dft_solver.grid_data_standard.quadrature_nodes
            quadrature_weights = atomic_dft_solver.grid_data_standard.quadrature_weights

            # Save the v_x_local and v_c_local on the uniform grid
            quadrature_nodes_uniform = final_result['uniform_grid']
            v_x_local_uniform = final_result['v_x_local_on_uniform_grid']
            v_c_local_uniform = final_result['v_c_local_on_uniform_grid']
            e_x_local_uniform = final_result['e_x_local_on_uniform_grid']
            e_c_local_uniform = final_result['e_c_local_on_uniform_grid']

            # Create output folder if it doesn't exist
            os.makedirs(write_folder_path, exist_ok=True)
            
            # Save data as txt files (original grid)
            np.savetxt(os.path.join(write_folder_path, "quadrature_nodes.txt"), quadrature_nodes)
            np.savetxt(os.path.join(write_folder_path, "quadrature_weights.txt"), quadrature_weights)
            np.savetxt(os.path.join(write_folder_path, "v_x.txt"), v_x_local)
            np.savetxt(os.path.join(write_folder_path, "v_c.txt"), v_c_local)
            
            # Save data as txt files (uniform grid)
            np.savetxt(os.path.join(write_folder_path, "quadrature_nodes_uniform.txt"), quadrature_nodes_uniform)
            np.savetxt(os.path.join(write_folder_path, "v_x_uniform.txt"), v_x_local_uniform)
            np.savetxt(os.path.join(write_folder_path, "v_c_uniform.txt"), v_c_local_uniform)

            # Save full_eigen_energies / full_l_terms (save_full_spectrum) and full_orbitals (save_full_orbitals)
            if save_full_spectrum:
                if final_result.get('full_eigen_energies') is not None:
                    np.savetxt(os.path.join(write_folder_path, "full_eigen_energies.txt"), final_result['full_eigen_energies'])
                if final_result.get('full_l_terms') is not None:
                    np.savetxt(os.path.join(write_folder_path, "full_l_terms.txt"), final_result['full_l_terms'])
            if save_full_orbitals and final_result.get('full_orbitals') is not None:
                np.savetxt(os.path.join(write_folder_path, "full_orbitals.txt"), final_result['full_orbitals'])

            if compute_energy_density:
                # Save local XC energy density
                np.savetxt(os.path.join(write_folder_path, "e_x.txt"), e_x_local)
                np.savetxt(os.path.join(write_folder_path, "e_c.txt"), e_c_local)
                # Save local XC energy density on uniform grid
                np.savetxt(os.path.join(write_folder_path, "e_x_uniform.txt"), e_x_local_uniform)
                np.savetxt(os.path.join(write_folder_path, "e_c_uniform.txt"), e_c_local_uniform)

            # Save output log for this folder (only if requested)
            # For nested calls, we save a snapshot of the current buffer content
            # The main function will save the complete log at the end
            # Intermediate subfolders should not save out.txt to avoid clutter
            if save_output_log:
                os.makedirs(write_folder_path, exist_ok=True)
                current_output = output_buffer.getvalue()
                with open(os.path.join(write_folder_path, "out.txt"), "w") as f:
                    f.write(current_output)

            # Restore stdout if needed
            if should_restore_stdout:
                sys.stdout = original_stdout
                    
        except Exception as e:
            # Save output log even on error (only if requested)
            # Intermediate subfolders should not save out.txt to avoid clutter
            if save_output_log:
                os.makedirs(write_folder_path, exist_ok=True)
                current_output = output_buffer.getvalue()
                with open(os.path.join(write_folder_path, "out.txt"), "w") as f:
                    f.write(current_output)
                    f.write(f"\n\nError occurred: {str(e)}\n")    
            
            # Restore stdout if needed
            if should_restore_stdout:
                sys.stdout = original_stdout
            raise
    
        return (v_x_local, v_c_local, e_x_local, e_c_local), atomic_dft_solver


    @staticmethod
    def _get_configuration_folder_path(
        directory_path : str,
        index          : int,
    ) -> str:
        """
        Get configuration folder path with backward compatibility.
        
        Tries to find existing folder in this order:
        1. atom_XXX (old format, for backward compatibility)
        2. configuration_XXX (new format)
        If neither exists, returns new format path (for creating new folders).
        
        Args:
            directory_path : Root directory of the dataset
            index          : Configuration index (configuration_index or atomic_number)
        
        Returns:
            str: Path to configuration folder
        """
        folder_old = os.path.join(directory_path, f"atom_{index:03d}")
        folder_new = os.path.join(directory_path, f"configuration_{index:03d}")
        
        if os.path.exists(folder_old):
            return folder_old
        elif os.path.exists(folder_new):
            return folder_new
        else:
            # Default to new format for new folders
            return folder_new

    @classmethod
    def forward_pass_single_atom_data(
        cls,
        # Required arguments
        atomic_number              : int | float,
        n_electrons                : int | float,
        scf_xc_functional          : str,
        forward_pass_xc_functional : str,
        directory_path             : str,
        configuration_index        : Optional[int],

        # Arguments controlling the contents of the dataset
        compute_energy_density     : bool,
        process_intermediate       : bool,
        save_full_spectrum         : bool,
        save_full_orbitals         : bool,

        # Arguments controlling the generation process
        # Grid, basis, and mesh parameters
        domain_size                : float,
        finite_elements_number     : int,
        polynomial_order           : int,
        quadrature_point_number    : int,
        oep_basis_number           : int,
        mesh_type                  : str,
        mesh_concentration         : float,
        mesh_spacing               : float,

        # Debugging and verbose parameters
        verbose                    : bool,

    ) -> XCDataType:
    
        """
        Perform forward pass for a single folder.
        
        Required arguments
        ------------------
        `atomic_number` : int | float
            Atomic number.
        `n_electrons` : int | float
            Number of electrons (required, must not be None).
        `scf_xc_functional` : str
            SCF XC functional (used for full SCF calculation to convergence).
        `forward_pass_xc_functional` : str
            Forward pass XC functional (performs forward pass based on SCF results).
        `directory_path` : str
            Directory path.
        `configuration_index` : int | None
            Index of the atom in the list. If None, uses atomic_number.

        Arguments controlling the contents of the dataset
        -------------------------------------------------
        `compute_energy_density` : bool
            Whether to compute energy density.
        `process_intermediate` : bool
            Whether to process intermediate information.
        `save_full_spectrum` : bool
            If True, save full_eigen_energies and full_l_terms (default: False).
        `save_full_orbitals` : bool
            If True, save full_orbitals (default: False, to save storage).

        Arguments controlling the generation process
        ---------------------------------------------
        Grid, basis, and mesh parameters
        `domain_size` : float
            Radial computational domain size in atomic units.
        `finite_elements_number` : int
            Number of finite elements.
        `polynomial_order` : int
            Polynomial order of basis functions.
        `quadrature_point_number` : int
            Number of quadrature points.
        `oep_basis_number` : int
            OEP basis number.
        `mesh_type` : str
            Mesh distribution type.
        `mesh_concentration` : float
            Mesh concentration parameter.
        `mesh_spacing` : float
            Output uniform mesh spacing.

        Debugging and verbose parameters
        --------------------------------
        `verbose` : bool
            Whether to print information during execution.

        Returns
        -------
        `v_x_local` : np.ndarray
            Local XC potential.
        `v_c_local` : np.ndarray
            Local XC potential.
        `e_x_local` : np.ndarray | None
            Local XC energy density.
        `e_c_local` : np.ndarray | None
            Local XC energy density.
        """

        # Type checks only, other checks are done in the AtomicDFTSolver class
        if not isinstance(atomic_number, (int, float)):
            raise TypeError(ATOMIC_NUMBER_NOT_INT_OR_FLOAT_ERROR.format(atomic_number))
        if n_electrons is None:
            raise ValueError(N_ELECTRONS_IS_NONE_ERROR)
        if not isinstance(n_electrons, (int, float)):
            raise TypeError(N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(n_electrons))
        if not isinstance(domain_size, float):
            raise TypeError(DOMAIN_SIZE_NOT_FLOAT_ERROR.format(domain_size))
        if not isinstance(finite_elements_number, int):
            raise TypeError(FINITE_ELEMENTS_NUMBER_NOT_INTEGER_ERROR.format(finite_elements_number))
        if not isinstance(polynomial_order, int):
            raise TypeError(POLYNOMIAL_ORDER_NOT_INTEGER_ERROR.format(polynomial_order))
        if not isinstance(quadrature_point_number, int):
            raise TypeError(QUADRATURE_POINT_NUMBER_NOT_INTEGER_ERROR.format(quadrature_point_number))
        if not isinstance(oep_basis_number, int):
            raise TypeError(OEP_BASIS_NUMBER_NOT_INTEGER_ERROR.format(oep_basis_number))
        if not isinstance(verbose, bool):
            raise TypeError(VERBOSE_NOT_BOOL_ERROR.format(verbose))
        if not isinstance(mesh_type, str):
            raise TypeError(MESH_TYPE_NOT_STRING_ERROR.format(mesh_type))
        if not isinstance(mesh_concentration, float):
            raise TypeError(MESH_CONCENTRATION_NOT_FLOAT_ERROR.format(mesh_concentration))
        if not isinstance(scf_xc_functional, str):
            raise TypeError(SCF_XC_FUNCTIONAL_NOT_STRING_ERROR.format(scf_xc_functional))
        if not isinstance(forward_pass_xc_functional, str):
            raise TypeError(FORWARD_PASS_XC_FUNCTIONAL_NOT_STRING_ERROR.format(forward_pass_xc_functional))
        if not isinstance(directory_path, str):
            raise TypeError(DIRECTORY_PATH_NOT_STRING_ERROR.format(directory_path))
        if not isinstance(process_intermediate, bool):
            raise TypeError(PROCESS_INTERMEDIATE_NOT_BOOL_ERROR.format(process_intermediate))
        if not isinstance(save_full_spectrum, bool):
            raise TypeError(SAVE_FULL_SPECTRUM_NOT_BOOL_ERROR.format(save_full_spectrum))
        if not isinstance(save_full_orbitals, bool):
            raise TypeError(SAVE_FULL_ORBITALS_NOT_BOOL_ERROR.format(save_full_orbitals))
        
        # Determine configuration index and get folder path
        if configuration_index is None:
            if not float(atomic_number).is_integer():
                raise ValueError(ATOMIC_NUMBER_NOT_INTEGER_WHEN_USED_AS_INDEX_ERROR.format(atomic_number))
            configuration_index_value = int(atomic_number)
        else:
            configuration_index_value = configuration_index
        
        config_folder = cls._get_configuration_folder_path(directory_path, configuration_index_value)

        read_folder_path = os.path.join(config_folder, scf_xc_functional.lower())
        write_folder_path = os.path.join(config_folder, forward_pass_xc_functional.lower())

        # Setup output redirection to both terminal and buffer
        # When verbose=False, only capture to buffer (for out.txt), not to terminal
        output_buffer = StringIO()
        original_stdout = sys.stdout
        tee_output = TeeOutput(original_stdout, output_buffer, verbose=verbose)
        sys.stdout = tee_output

        try:
            # Perform forward pass for main folder
            print(f"Performing forward pass for main folder: {read_folder_path}")
            # Always use verbose=True in _forward_pass_single_folder to ensure detailed output is captured
            # TeeOutput will control whether output goes to terminal based on verbose parameter
            xc_data, atomic_dft_solver = \
                cls._forward_pass_single_folder(
                    # Required arguments
                    atomic_number           = atomic_number,
                    n_electrons             = n_electrons,
                    xc_functional           = forward_pass_xc_functional,
                    read_folder_path        = read_folder_path,
                    write_folder_path       = write_folder_path,

                    # Arguments controlling the contents of the dataset
                    compute_energy_density  = compute_energy_density,
                    save_full_spectrum      = save_full_spectrum,
                    save_full_orbitals      = save_full_orbitals,

                    # Arguments controlling the generation process
                    # Grid, basis, and mesh parameters
                    domain_size             = domain_size,
                    finite_elements_number  = finite_elements_number,
                    polynomial_order        = polynomial_order,
                    quadrature_point_number = quadrature_point_number,
                    oep_basis_number        = oep_basis_number,
                    mesh_type               = mesh_type,
                    mesh_concentration      = mesh_concentration,
                    mesh_spacing            = mesh_spacing,

                    # Debugging and verbose parameters
                    # Always use verbose=True to capture detailed output in buffer
                    # TeeOutput will control whether output goes to terminal
                    verbose                 = True,  # Always True to capture detailed output in buffer

                    # Optional parameters control
                    output_buffer           = output_buffer,
                    original_stdout         = original_stdout,
                    atomic_dft_solver       = None,
                )
            v_x_local, v_c_local, e_x_local, e_c_local = xc_data

            # check for outer iteration subfolders
            if process_intermediate and os.path.exists(read_folder_path):
                outer_iter_folders = [f for f in os.listdir(read_folder_path) 
                                      if os.path.isdir(os.path.join(read_folder_path, f)) 
                                      and f.startswith("outer_iter_")]
                outer_iter_folders.sort()
                
                if len(outer_iter_folders) > 0:
                    print(f"Found {len(outer_iter_folders)} outer iteration folders, processing...")
                    
                    # process each outer iteration folder
                    for outer_iter_folder in outer_iter_folders:
                        outer_iter_read_path  = os.path.join(read_folder_path, outer_iter_folder)
                        outer_iter_write_path = os.path.join(write_folder_path, outer_iter_folder)
                        
                        print(f"Processing {outer_iter_folder}...")
                        try:
                            # Reuse the same solver instance to avoid creating new one (reduces verbose output)
                            cls._forward_pass_single_folder(
                                # Required arguments
                                atomic_number           = atomic_number,
                                n_electrons             = n_electrons,
                                xc_functional           = forward_pass_xc_functional,
                                read_folder_path        = outer_iter_read_path,  # read from intermediate sub-folder
                                write_folder_path       = outer_iter_write_path, # write to intermediate sub-folder

                                # Arguments controlling the contents of the dataset
                                compute_energy_density  = compute_energy_density,
                                save_full_spectrum      = save_full_spectrum,
                                save_full_orbitals      = save_full_orbitals,

                                # Arguments controlling the generation process
                                # Grid, basis, and mesh parameters
                                domain_size             = domain_size,
                                finite_elements_number  = finite_elements_number,
                                polynomial_order        = polynomial_order,
                                quadrature_point_number = quadrature_point_number,
                                oep_basis_number        = oep_basis_number,
                                mesh_type               = mesh_type,
                                mesh_concentration      = mesh_concentration,
                                mesh_spacing            = mesh_spacing,

                                # Debugging and verbose parameters
                                verbose                 = False,                 # Disable information printing for intermediate folders to reduce output

                                # Optional parameters control
                                output_buffer           = output_buffer,
                                original_stdout         = original_stdout,
                                atomic_dft_solver       = atomic_dft_solver,     # Reuse solver
                                save_output_log         = False                  # Don't save out.txt in intermediate subfolders
                            )
                            print(f"  ✓ Successfully processed {outer_iter_folder} \n")
                        except Exception as e:
                            print(f"  ✗ Error processing {outer_iter_folder}: {str(e)} \n")
                            # Continue with other folders even if one fails
                            continue
            
            # Restore stdout and save final log
            sys.stdout = original_stdout
            os.makedirs(write_folder_path, exist_ok=True)
            with open(os.path.join(write_folder_path, "out.txt"), "w") as f:
                f.write(output_buffer.getvalue())

            print(f"Successfully performed forward pass for atom {atomic_number} in folder '{write_folder_path}'")

        except Exception as e:
            sys.stdout = original_stdout
            error_summary, error_traceback = format_error_message(e, f"Error performing forward pass for atom {atomic_number}")
            
            print(error_summary)
            print(f"Traceback:\n{error_traceback}")

            os.makedirs(write_folder_path, exist_ok=True)
            with open(os.path.join(write_folder_path, "out.txt"), "w") as f:
                f.write(output_buffer.getvalue())
                f.write(f"\n\n{error_summary}\n")
                f.write(f"\nTraceback:\n{error_traceback}\n")
            
            raise

        return v_x_local, v_c_local, e_x_local, e_c_local


    @classmethod
    def generate_data(cls, 
        # Required arguments
        data_root                   : str,
        atomic_number_list          : list[int | float], 
        n_electrons_list            : Optional[list[int | float]],
        use_oep                     : bool,
        scf_xc_functional           : str,
        forward_pass_xc_functionals : Optional[str | list[str]],

        # Arguments controlling the contents of the dataset
        save_energy_density         : bool,
        save_intermediate           : bool,
        save_full_spectrum          : bool,
        save_full_orbitals          : bool,
        save_derivative_matrix      : bool,
        start_configuration_index   : int,

        # Arguments controlling the generation process
        # Grid, basis, and mesh parameters
        domain_size                 : float,
        finite_elements_number      : int,
        polynomial_order            : int,
        quadrature_point_number     : int,
        oep_basis_number            : int,
        mesh_type                   : str,
        mesh_concentration          : float,
        mesh_spacing                : float,

        # SCF convergence parameters
        scf_tolerance               : float,
        max_scf_iterations          : int,
        max_scf_iterations_outer    : Optional[int],
        use_pulay_mixing            : bool,
        use_preconditioner          : bool,
        pulay_mixing_parameter      : float,
        pulay_mixing_history        : int,
        pulay_mixing_frequency      : int,
        linear_mixing_alpha1        : float,
        linear_mixing_alpha2        : float,

        # Advanced functional parameters
        hybrid_mixing_parameter           : float,
        frequency_quadrature_point_number : int,
        angular_momentum_cutoff           : int,
        double_hybrid_flag                : bool,
        oep_mixing_parameter              : float,
        enable_parallelization            : bool,

        # Debugging and verbose parameters
        verbose                     : bool,
    ):
        """
        Generate atomic dataset, based on the source and target XC functionals.
        
        Required arguments
        ------------------
        `data_root` : str
            Root directory of the dataset.
        `atomic_number_list` : list[int | float]
            List of atomic numbers to generate data for.
        `n_electrons_list` : list[int | float] | None
            List of number of electrons. If None, defaults to atomic_number_list.
        `use_oep` : bool
            Whether to use OEP for the SCF calculation.
        `scf_xc_functional` : str
            SCF XC functional (used for full SCF calculation to convergence).
        `forward_pass_xc_functionals` : str | list[str] | None
            Forward pass XC functional(s). If not None, will perform forward pass for each functional based on SCF results.

        Dataset content control
        -----------------------
        `save_energy_density` : bool
            Whether to save energy density.
        `save_intermediate` : bool
            Whether to save intermediate information.
        `save_full_spectrum` : bool
            Whether to save full_eigen_energies and full_l_terms.
        `save_full_orbitals` : bool
            Whether to save full_orbitals (default: False, to save storage).
        `save_derivative_matrix` : bool
            Whether to save derivative matrix. Most systems have the same derivative matrix when using
            the same grid/basis/mesh parameters, so a shared derivative matrix is saved at the dataset root.
            If an atom's derivative matrix differs from the shared one, it is saved locally and recorded in meta.json.
        `start_configuration_index` : int
            Starting configuration index for generated folders. Configuration folders will be named
            configuration_XXX starting from this index.
            For example, if start_configuration_index=5, the first atom will be saved as configuration_005,
            the second as configuration_006, etc.

        Grid, basis, and mesh parameters
        --------------------------------
        `domain_size` : float
            Radial computational domain size in atomic units.
        `finite_elements_number` : int
            Number of finite elements in the computational domain.
        `polynomial_order` : int
            Polynomial order of basis functions.
        `quadrature_point_number` : int
            Number of quadrature points for numerical integration.
        `oep_basis_number` : int
            Basis size used in OEP calculations.
        `mesh_type` : str
            Mesh distribution type ('exponential', 'polynomial', 'uniform').
        `mesh_concentration` : float
            Mesh concentration parameter.
        `mesh_spacing` : float
            Output uniform mesh spacing.

        SCF convergence parameters
        --------------------------
        `scf_tolerance` : float
            SCF convergence tolerance.
        `max_scf_iterations` : int | None
            Maximum number of inner SCF iterations. If None, uses default (500).
        `max_scf_iterations_outer` : int | None
            Maximum number of outer SCF iterations (for functionals requiring outer loop like HF, EXX, RPA, PBE0).
            If None, uses default (50 when needed, otherwise not used).
        `use_pulay_mixing` : bool
            True for Pulay mixing, False for linear mixing.
        `use_preconditioner` : bool
            Whether to use preconditioner for SCF convergence.
        `pulay_mixing_parameter` : float
            Pulay mixing parameter.
        `pulay_mixing_history` : int
            Pulay mixing history.
        `pulay_mixing_frequency` : int
            Pulay mixing frequency.
        `linear_mixing_alpha1` : float
            Linear mixing parameter (alpha_1).
        `linear_mixing_alpha2` : float
            Linear mixing parameter (alpha_2).

        Advanced functional parameters
        ------------------------------
        `hybrid_mixing_parameter` : float
            Mixing parameter for hybrid/double-hybrid functionals.
        `frequency_quadrature_point_number` : int
            Number of frequency quadrature points for RPA calculations.
        `angular_momentum_cutoff` : int
            Maximum angular momentum quantum number to include.
        `double_hybrid_flag` : bool
            Flag for double-hybrid functional methods.
        `oep_mixing_parameter` : float
            Scaling parameter (λ) for OEP exchange/correlation potentials.
        `enable_parallelization` : bool
            Flag for parallelization of RPA calculations.

        Debugging and verbose parameters
        --------------------------------
        `verbose` : bool
            Whether to print information during execution.
        """
        # atomic_number_list
        if not isinstance(atomic_number_list, list):
            raise TypeError(ATOMIC_NUMBER_LIST_NOT_LIST_OF_INTEGERS_ERROR.format(atomic_number_list))
        for atomic_number in atomic_number_list:
            if not isinstance(atomic_number, (int, float)):
                raise TypeError(ATOMIC_NUMBER_NOT_INT_OR_FLOAT_ERROR.format(atomic_number))
    
        # scf_xc_functional
        if not isinstance(scf_xc_functional, str):
            raise TypeError(SCF_XC_FUNCTIONAL_NOT_STRING_ERROR.format(scf_xc_functional))
        
        # start_configuration_index
        if not isinstance(start_configuration_index, int):
            raise TypeError(START_CONFIGURATION_INDEX_NOT_INTEGER_ERROR.format(start_configuration_index))
        if start_configuration_index < 1:
            raise ValueError(START_CONFIGURATION_INDEX_NOT_POSITIVE_ERROR.format(start_configuration_index))
        
        # forward_pass_xc_functionals
        if forward_pass_xc_functionals is None:
            forward_pass_xc_functional_list = []
        elif isinstance(forward_pass_xc_functionals, str):
            forward_pass_xc_functional_list = [forward_pass_xc_functionals,]
        elif isinstance(forward_pass_xc_functionals, list):
            forward_pass_xc_functional_list = forward_pass_xc_functionals
        else:
            raise TypeError(FORWARD_PASS_XC_FUNCTIONAL_NOT_STRING_OR_LIST_OF_STRINGS_ERROR.format(forward_pass_xc_functionals))

        # n_electrons_list
        if n_electrons_list is not None:
            if not isinstance(n_electrons_list, list):
                raise TypeError(N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(n_electrons_list))
            if len(n_electrons_list) != len(atomic_number_list):
                raise ValueError(N_ELECTRONS_LIST_LENGTH_MISMATCH_ERROR.format(
                    len(n_electrons_list), len(atomic_number_list)
                ))
            for n_electrons in n_electrons_list:
                if not isinstance(n_electrons, (int, float)):
                    raise TypeError(N_ELECTRONS_NOT_INT_OR_FLOAT_ERROR.format(n_electrons))
        else:
            # default to charge neutral
            n_electrons_list = atomic_number_list 

        # Initialize shared derivative matrix for the dataset (if save_derivative_matrix is True)
        shared_derivative_matrix_path = None
        if save_derivative_matrix:
            shared_derivative_matrix_path = os.path.join(data_root, "derivative_matrix.npy")

        # Generate data for specified atomic numbers
        for idx, atomic_number in enumerate(atomic_number_list, 1):
            n_electrons = n_electrons_list[idx - 1]
            # Calculate configuration index: start from start_configuration_index
            configuration_index = start_configuration_index + idx - 1
            print("\n")
            print("=" * 70)
            print(f"Processing atom {atomic_number} ({idx}/{len(atomic_number_list)}, configuration_index={configuration_index})")
            print("=" * 70)

            try:
                cls.generate_single_atom_data(
                    # Required arguments
                    atomic_number            = atomic_number,
                    n_electrons              = n_electrons,
                    use_oep                  = use_oep,
                    directory_path           = data_root,
                    configuration_index      = configuration_index,

                    # Arguments controlling the contents of the dataset
                    save_energy_density      = save_energy_density,
                    save_intermediate        = save_intermediate,
                    save_full_spectrum       = save_full_spectrum,
                    save_full_orbitals       = save_full_orbitals,
                    save_derivative_matrix   = save_derivative_matrix,

                    # Arguments controlling the generation process
                    # Grid, basis, and mesh parameters
                    domain_size              = domain_size,
                    finite_elements_number   = finite_elements_number,
                    polynomial_order         = polynomial_order,
                    quadrature_point_number  = quadrature_point_number,
                    oep_basis_number         = oep_basis_number,
                    mesh_type                = mesh_type,
                    mesh_concentration       = mesh_concentration,
                    mesh_spacing             = mesh_spacing,

                    # SCF convergence parameters
                    xc_functional            = scf_xc_functional,
                    scf_tolerance            = scf_tolerance,
                    max_scf_iterations       = max_scf_iterations,
                    max_scf_iterations_outer = max_scf_iterations_outer,
                    use_pulay_mixing         = use_pulay_mixing,
                    use_preconditioner       = use_preconditioner,
                    pulay_mixing_parameter   = pulay_mixing_parameter,
                    pulay_mixing_history     = pulay_mixing_history,
                    pulay_mixing_frequency   = pulay_mixing_frequency,
                    linear_mixing_alpha1     = linear_mixing_alpha1,
                    linear_mixing_alpha2     = linear_mixing_alpha2,

                    # Advanced functional parameters
                    hybrid_mixing_parameter           = hybrid_mixing_parameter,
                    frequency_quadrature_point_number = frequency_quadrature_point_number,
                    angular_momentum_cutoff           = angular_momentum_cutoff,
                    double_hybrid_flag                = double_hybrid_flag,
                    oep_mixing_parameter              = oep_mixing_parameter,
                    enable_parallelization            = enable_parallelization,

                    # Debugging and verbose parameters
                    verbose = verbose,

                    # Derivative matrix handling
                    shared_derivative_matrix_path = shared_derivative_matrix_path,
                )
                for forward_pass_xc_functional in forward_pass_xc_functional_list:
                    cls.forward_pass_single_atom_data(
                        # Required arguments
                        atomic_number              = atomic_number,
                        n_electrons                = n_electrons,
                        scf_xc_functional          = scf_xc_functional,
                        forward_pass_xc_functional = forward_pass_xc_functional,
                        directory_path             = data_root,
                        configuration_index        = configuration_index,

                        # Arguments controlling the contents of the dataset
                        compute_energy_density     = save_energy_density,
                        process_intermediate       = save_intermediate,
                        save_full_spectrum         = save_full_spectrum,
                        save_full_orbitals         = save_full_orbitals,

                        # Arguments controlling the generation process
                        # Grid, basis, and mesh parameters
                        domain_size                = domain_size,
                        finite_elements_number     = finite_elements_number,
                        polynomial_order           = polynomial_order,
                        quadrature_point_number    = quadrature_point_number,
                        oep_basis_number           = oep_basis_number,
                        mesh_type                  = mesh_type,
                        mesh_concentration         = mesh_concentration,
                        mesh_spacing               = mesh_spacing,

                        # Debugging and verbose parameters
                        # Note: Advanced functional parameters are not used in forward pass
                        verbose                    = verbose,
                    )
            except Exception as e:
                error_summary, error_traceback = format_error_message(e, f"Error generating data for atom {atomic_number}")
                print(error_summary)
                print(f"Traceback:\n{error_traceback}")

                if (idx + 1) < len(atomic_number_list):
                    print("Continuing to next atom...")
                continue
        print()
        print("="*75)
        print("Data generation complete".center(75))
        print("="*75)

