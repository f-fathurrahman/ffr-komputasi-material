
from __future__ import annotations
import numpy as np
from typing import Optional

# Error Messages for initial parameters
USE_PULAY_MIXING_NOT_BOOL_ERROR = \
    "parameter 'use_pulay_mixing' must be a boolean, get {} instead."
USE_PRECONDITIONER_NOT_BOOL_ERROR = \
    "parameter 'use_preconditioner' must be a boolean, get {} instead."
USE_RHO_CLAMP_NOT_BOOL_ERROR = \
    "parameter 'use_rho_clamp' must be a boolean, get {} instead."
PULAY_MIXING_PARAMETER_NOT_FLOAT_ERROR = \
    "parameter 'pulay_mixing_parameter' must be a float, get {} instead."
PULAY_MIXING_PARAMETER_NOT_IN_ZERO_ONE_ERROR = \
    "parameter 'pulay_mixing_parameter' must be in [0, 1], get {} instead."
PULAY_MIXING_HISTORY_NOT_INTEGER_ERROR = \
    "parameter 'pulay_mixing_history' must be an integer, get {} instead."
PULAY_MIXING_HISTORY_NOT_GREATER_THAN_0_ERROR = \
    "parameter 'pulay_mixing_history' must be greater than 0, get {} instead."
PULAY_MIXING_FREQUENCY_NOT_INTEGER_ERROR = \
    "parameter 'pulay_mixing_frequency' must be an integer, get {} instead."
PULAY_MIXING_FREQUENCY_NOT_GREATER_THAN_0_ERROR = \
    "parameter 'pulay_mixing_frequency' must be greater than 0, get {} instead."
LINEAR_MIXING_ALPHA1_NOT_FLOAT_ERROR = \
    "parameter 'linear_mixing_alpha1' must be a float, get {} instead."
LINEAR_MIXING_ALPHA1_NOT_IN_ZERO_ONE_ERROR = \
    "parameter 'linear_mixing_alpha1' must be in [0, 1], get {} instead."
LINEAR_MIXING_ALPHA2_NOT_FLOAT_ERROR = \
    "parameter 'linear_mixing_alpha2' must be a float, get {} instead."
LINEAR_MIXING_ALPHA2_NOT_IN_ZERO_ONE_ERROR = \
    "parameter 'linear_mixing_alpha2' must be in [0, 1], get {} instead."
RHO_CLAMP_MINIMUM_NOT_FLOAT_ERROR = \
    "parameter 'rho_clamp_minimum' must be a float, get {} instead."
RHO_CLAMP_MINIMUM_NOT_GREATER_THAN_0_ERROR = \
    "parameter 'rho_clamp_minimum' must be greater than 0, get {} instead."

PRECONDITIONER_IS_NONE_WHEN_USE_PRECONDITIONER_IS_TRUE_ERROR = \
    "parameter 'preconditioner' must be not None when 'use_preconditioner' is True."


# Warning Messages for initial parameters
PULAY_MIXING_PARAMETER_NOT_NONE_WHEN_USE_PULAY_MIXING_IS_FALSE_WARNING = \
    "WARNING: parameter 'pulay_mixing_parameter' is not None when 'use_pulay_mixing' is False, so it will be ignored."
PULAY_MIXING_HISTORY_NOT_NONE_WHEN_USE_PULAY_MIXING_IS_FALSE_WARNING = \
    "WARNING: parameter 'pulay_mixing_history' is not None when 'use_pulay_mixing' is False, so it will be ignored."
PULAY_MIXING_FREQUENCY_NOT_NONE_WHEN_USE_PULAY_MIXING_IS_FALSE_WARNING = \
    "WARNING: parameter 'pulay_mixing_frequency' is not None when 'use_pulay_mixing' is False, so it will be ignored."
PRECONDITIONER_IS_NOT_NONE_WHEN_USE_PRECONDITIONER_IS_FALSE_WARNING = \
    "WARNING: parameter 'preconditioner' is not None when 'use_preconditioner' is False, so it will be ignored."



class Mixer:
    """
    Density mixing for SCF convergence.
    

    Supports:
    - Simple linear mixing: rho_new = rho_out (pulay_mixing=0)
    - Pulay mixing: DIIS-like method using history of residuals
    - Alternating linear mixing: odd/even steps with different coefficients


    References
    ----------
    Periodic Pulay method for robust and efficient convergence acceleration of self-consistent field iterations,
    Amartya S. Banerjee, Phanish Suryanarayana, John E. Pask,
    Chemical Physics Letters, Volume 647, 2016, Pages 31-35, ISSN 0009-2614, https://doi.org/10.1016/j.cplett.2016.01.033.


    Parameters
    ----------
    use_pulay_mixing : bool
        Whether to use Pulay mixing, if False, use linear mixing.
    use_preconditioner : bool
        Whether to use preconditioner. Can be used with both Pulay mixing and linear mixing.
    pulay_mixing_parameter : float
        Pulay mixing parameter
    pulay_mixing_history : int
        Number of previous iterations to keep for Pulay mixing
    pulay_mixing_frequency : int
        Apply Pulay mixing every 'frequency' iterations
    linear_mixing_alpha1 : float
        Linear mixing parameter (alpha21) for odd/even steps
    linear_mixing_alpha2 : float
        Linear mixing parameter (alpha22) for even/odd steps
    rho_clamp_minimum : float
        Minimum density value to avoid negative densities


    Algorithm Walkthrough:
    -----------------------
    history = 4, frequency = 2

    Initialization:
        rho_in  = [*, o, o, o, o]
        rho_out = [o, o, o, o, o]

    Step 1: run = 0
        Solve KS equation:
            rho_in  = [x, o, o, o, o]
            rho_out = [*, o, o, o, o]
        Mix: (linear mixing, history not full)
            rho_in  = [x, *, o, o, o]
            rho_out = [x, o, o, o, o]
        
    Step 2: run = 1
        Solve KS equation:
            rho_in  = [x, x, o, o, o]
            rho_out = [x, *, o, o, o]
        Mix: (Pulay mixing, (run+1) % frequency == 0)
            rho_in  = [x, x, *, o, o]
            rho_out = [x, x, o, o, o]
        
    Step 3: run = 2
        Solve KS equation:
            rho_in  = [x, x, x, o, o]
            rho_out = [x, x, *, o, o]
        Mix: (linear mixing, (run+1) % frequency != 0)
            rho_in  = [x, x, x, *, o]
            rho_out = [x, x, x, o, o]
    
    Step 4: run = 3
        Solve KS equation:
            rho_in  = [x, x, x, x, o]
            rho_out = [x, x, x, *, o]
        Mix: (Pulay mixing, (run+1) % frequency == 0, history full)
            rho_in  = [x, x, x, x, *]
            rho_out = [x, x, x, x, o]

    Step 5: run = 4
        Solve KS equation:
            rho_in  = [x, x, x, x, x]
            rho_out = [x, x, x, x, *]
        Mix: (linear mixing, (run+1) % frequency != 0)
            rho_in  = x, [x, x, x, x, *]
            rho_out = x, [x, x, x, x, o]
    
    Step 6: run = 5
        Solve KS equation:
            rho_in  = [x, x, x, x, x]
            rho_out = [x, x, x, x, *]
        Mix: (Pulay mixing, (run+1) % frequency == 0)
            rho_in  = x, [x, x, x, x, *]
            rho_out = x, [x, x, x, x, o]
    
    Step 7: run = 6
        Same as Step 5 and 6, til convergence.
    """


    def __init__(
        self, 
        use_pulay_mixing       : bool  = True,
        use_preconditioner     : bool  = True,
        pulay_mixing_parameter : Optional[float] = None, # default is 1.0 if use_preconditioner else 0.45
        pulay_mixing_history   : Optional[int]   = None, # default is 7 if use_preconditioner else 11
        pulay_mixing_frequency : Optional[int]   = None, # default is 3 if use_preconditioner else 1
        linear_mixing_alpha1   : Optional[float] = None, # default is 0.75 if use_pulay_mixing else 0.7
        linear_mixing_alpha2   : Optional[float] = None, # default is 0.95 if use_pulay_mixing else 1.0
        rho_clamp_minimum      : Optional[float] = None, # default is 1e-20
    ):
        self.use_pulay_mixing       = use_pulay_mixing
        self.use_preconditioner     = use_preconditioner
        self.pulay_mixing_parameter = pulay_mixing_parameter 
        self.pulay_mixing_history   = pulay_mixing_history
        self.pulay_mixing_frequency = pulay_mixing_frequency
        self.linear_mixing_alpha1   = linear_mixing_alpha1  # Odd step mixing
        self.linear_mixing_alpha2   = linear_mixing_alpha2  # Even step mixing
        self.use_rho_clamp          = True # default is True, use self.set_use_rho_clamp to change
        self.rho_clamp_minimum      = rho_clamp_minimum
        self._set_and_check_initial_parameters()
        
        # History storage: (N_q, history+1)
        self.iteration_count : int = 0
        self.n_points        : Optional[int]        = None
        self.rho_in_store    : Optional[np.ndarray] = None
        self.rho_out_store   : Optional[np.ndarray] = None


    def _set_and_check_initial_parameters(self) -> None:
        """
        Check the initial parameters.
        """
        # use pulay mixing flag
        if self.use_pulay_mixing is None:
            self.use_pulay_mixing = True
        if self.use_pulay_mixing in [0, 1]:
            self.use_pulay_mixing = False if self.use_pulay_mixing == 0 else True
        assert isinstance(self.use_pulay_mixing, bool), \
            USE_PULAY_MIXING_NOT_BOOL_ERROR.format(type(self.use_pulay_mixing))


        # use preconditioner flag
        if self.use_preconditioner is None:
            self.use_preconditioner = True if self.use_pulay_mixing else False
        if self.use_preconditioner in [0, 1]:
            self.use_preconditioner = False if self.use_preconditioner == 0 else True
        assert isinstance(self.use_preconditioner, bool), \
            USE_PRECONDITIONER_NOT_BOOL_ERROR.format(type(self.use_preconditioner))

        # pulay mixing parameter
        if self.use_pulay_mixing:
            if self.pulay_mixing_parameter is None:
                self.pulay_mixing_parameter = 1.0 if self.use_preconditioner else 0.45
            try:
                self.pulay_mixing_parameter = float(self.pulay_mixing_parameter)
            except:
                raise ValueError(PULAY_MIXING_PARAMETER_NOT_FLOAT_ERROR.format(type(self.pulay_mixing_parameter)))
            assert isinstance(self.pulay_mixing_parameter, float), \
                PULAY_MIXING_PARAMETER_NOT_FLOAT_ERROR.format(type(self.pulay_mixing_parameter))
            assert 0.0 < self.pulay_mixing_parameter <= 1.0, \
                PULAY_MIXING_PARAMETER_NOT_IN_ZERO_ONE_ERROR.format(self.pulay_mixing_parameter)
        else:
            if self.pulay_mixing_parameter is not None:
                print(PULAY_MIXING_PARAMETER_NOT_NONE_WHEN_USE_PULAY_MIXING_IS_FALSE_WARNING)
                self.pulay_mixing_parameter = None


        # pulay mixing history
        if self.use_pulay_mixing:
            if self.pulay_mixing_history is None:
                self.pulay_mixing_history = 7 if self.use_preconditioner else 11
            assert isinstance(self.pulay_mixing_history, int), \
                PULAY_MIXING_HISTORY_NOT_INTEGER_ERROR.format(type(self.pulay_mixing_history))
            assert self.pulay_mixing_history > 0, \
                PULAY_MIXING_HISTORY_NOT_GREATER_THAN_0_ERROR.format(self.pulay_mixing_history)
        else:
            if self.pulay_mixing_history is not None:
                print(PULAY_MIXING_HISTORY_NOT_NONE_WHEN_USE_PULAY_MIXING_IS_FALSE_WARNING)
                self.pulay_mixing_history = None


        # pulay mixing frequency
        if self.use_pulay_mixing:
            if self.pulay_mixing_frequency is None:
                self.pulay_mixing_frequency = 3 if self.use_preconditioner else 1
            assert isinstance(self.pulay_mixing_frequency, int), \
                PULAY_MIXING_FREQUENCY_NOT_INTEGER_ERROR.format(type(self.pulay_mixing_frequency))
            assert self.pulay_mixing_frequency > 0, \
                PULAY_MIXING_FREQUENCY_NOT_GREATER_THAN_0_ERROR.format(self.pulay_mixing_frequency)
        else:   
            if self.pulay_mixing_frequency is not None:
                print(PULAY_MIXING_FREQUENCY_NOT_NONE_WHEN_USE_PULAY_MIXING_IS_FALSE_WARNING)
                self.pulay_mixing_frequency = None


        # linear mixing alpha1
        if self.linear_mixing_alpha1 is None:
            self.linear_mixing_alpha1 = 0.75 if self.use_pulay_mixing else 0.7
        assert isinstance(self.linear_mixing_alpha1, float), \
            LINEAR_MIXING_ALPHA1_NOT_FLOAT_ERROR.format(type(self.linear_mixing_alpha1))
        assert 0.0 <= self.linear_mixing_alpha1 <= 1.0, \
            LINEAR_MIXING_ALPHA1_NOT_IN_ZERO_ONE_ERROR.format(self.linear_mixing_alpha1)

        # linear mixing alpha2
        if self.linear_mixing_alpha2 is None:
            self.linear_mixing_alpha2 = 0.95 if self.use_pulay_mixing else 1.0
        assert isinstance(self.linear_mixing_alpha2, float), \
            LINEAR_MIXING_ALPHA2_NOT_FLOAT_ERROR.format(type(self.linear_mixing_alpha2))
        assert 0.0 <= self.linear_mixing_alpha2 <= 1.0, \
            LINEAR_MIXING_ALPHA2_NOT_IN_ZERO_ONE_ERROR.format(self.linear_mixing_alpha2)

        # rho clamp minimum
        if self.rho_clamp_minimum is None:
            self.rho_clamp_minimum = 1e-20
        assert isinstance(self.rho_clamp_minimum, float), \
            RHO_CLAMP_MINIMUM_NOT_FLOAT_ERROR.format(type(self.rho_clamp_minimum))
        assert self.rho_clamp_minimum > 0, \
            RHO_CLAMP_MINIMUM_NOT_GREATER_THAN_0_ERROR.format(self.rho_clamp_minimum)


    def reset(self) -> None:
        """
        Reset the mixer state.
        
        Clears all history buffers and resets iteration counter.
        Called at the start of each SCF calculation.
        """
        self.iteration_count = 0
        self.n_points        = None
        self.rho_in_store    = None
        self.rho_out_store   = None

    def set_use_rho_clamp(self, use_rho_clamp: bool = True) -> None:
        """
        Enable/disable rho clamping during mixing.

        Parameters
        ----------
        use_rho_clamp : bool, default=True
            If True, clamp rho <= 0 to rho_clamp_minimum.
        """
        if use_rho_clamp in [0, 1]:
            use_rho_clamp = False if use_rho_clamp == 0 else True
        assert isinstance(use_rho_clamp, bool), \
            USE_RHO_CLAMP_NOT_BOOL_ERROR.format(type(use_rho_clamp))
        self.use_rho_clamp = use_rho_clamp

    @staticmethod
    def _safe_solve_linear_system(A: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        Solve A x = b with a fallback for unsupported extended dtypes.
        """
        try:
            return np.linalg.solve(A, b)
        except TypeError:
            return np.linalg.solve(np.asarray(A, dtype=np.float64), np.asarray(b, dtype=np.float64))

    
    def mix(
        self, 
        rho_in         : np.ndarray, 
        rho_out        : np.ndarray,
        preconditioner : Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Mix input and output densities to produce next input density.
        
        Implements the exact logic from the original code, including:
        - History management with rolling window
        - Pulay mixing every 'frequency' iterations
        - Alternating linear mixing (odd/even steps)
        - Negative density correction
        
        Parameters
        ----------
        rho_in : np.ndarray
            Input density for current iteration, shape (n_points,)
        rho_out : np.ndarray
            Output density from KS solve, shape (n_points,)
        preconditioner : Optional[np.ndarray]
            Preconditioner matrix, shape (n_points, n_points), default is None
            If use_preconditioner is True, preconditioner is required.

        Returns
        -------
        rho_next : np.ndarray
            Mixed density for next iteration
        
        Example
        --------
        >>> mixer = Mixer()
        >>> rho = np.zeros(10)
        >>> while not converged:
        >>>     rho_out, preconditioner = KS_solve(rho)
        >>>     rho = mixer.mix(rho, rho_out, preconditioner)
        """

        if self.use_preconditioner:
            assert preconditioner is not None, \
                PRECONDITIONER_IS_NONE_WHEN_USE_PRECONDITIONER_IS_TRUE_ERROR
        else:
            if preconditioner is not None:
                print(PRECONDITIONER_IS_NOT_NONE_WHEN_USE_PRECONDITIONER_IS_FALSE_WARNING)
                preconditioner = None

        if self.use_pulay_mixing:
            # Initialize history storage on first call
            if self.rho_in_store is None:
                self.n_points      = rho_in.shape[0]
                self.rho_in_store  = np.zeros((self.n_points, self.pulay_mixing_history + 1))
                self.rho_out_store = np.zeros((self.n_points, self.pulay_mixing_history + 1))
                self.rho_in_store[:, 0] = rho_in
                self.iteration_count = 0
            
            runs = self.iteration_count
            
            # Store rho_out in appropriate column
            if runs <= self.pulay_mixing_history:
                self.rho_out_store[:, runs] = rho_out
            else:
                # self.rho_out_store = np.roll(self.rho_out_store, -1, axis=1)
                self.rho_out_store[:, -1] = rho_out
            
            # Pulay + alternating linear mixing
            rho_next = self._hybrid_mix(runs, preconditioner)
            
            # Store next rho_in for next iteration
            if runs < self.pulay_mixing_history:
                self.rho_in_store[:, runs + 1] = rho_next
            else:
                # Roll history and store in last column
                self.rho_in_store  = np.roll(self.rho_in_store,  -1, axis=1)
                self.rho_out_store = np.roll(self.rho_out_store, -1, axis=1)
                self.rho_in_store[:, -1] = rho_next
        else:
            # Linear mixing
            runs = self.iteration_count
            rho_next = self._linear_mix(rho_in, rho_out, runs, preconditioner)

        # Increment iteration count
        self.iteration_count += 1

        return rho_next
    

    def _linear_mix(
        self, 
        rho_in         : np.ndarray, 
        rho_out        : np.ndarray, 
        runs           : int, 
        preconditioner : Optional[np.ndarray]
    ) -> np.ndarray:
        """
        Linear mixing for early and stable phases.
        
        Parameters
        ----------
        rho_in : np.ndarray
            Input density for current iteration, shape (n_points,)
        rho_out : np.ndarray
            Output density from KS solve, shape (n_points,)
        runs : int
            Current iteration number (before increment in reference code)
        preconditioner : Optional[np.ndarray]
            Preconditioner matrix, shape (n_points, n_points), default is None
            Required when use_preconditioner is True.
        
        Returns
        -------
        rho_next : np.ndarray
            Mixed density        
        """

        if (runs + 1) % 2 != 0:  # equivalent to: runs % 2 == 0
            alpha = self.linear_mixing_alpha1
        else:
            alpha = self.linear_mixing_alpha2
        
        rho_residual_prev = rho_out - rho_in
        
        if self.use_preconditioner:
            assert preconditioner is not None, \
                PRECONDITIONER_IS_NONE_WHEN_USE_PRECONDITIONER_IS_TRUE_ERROR
            rho_residual_prev = np.linalg.solve(preconditioner, rho_residual_prev)
        
        return rho_in + alpha * rho_residual_prev



    def _hybrid_mix(self, runs: int, preconditioner: Optional[np.ndarray]) -> np.ndarray:
    
        """
        Hybrid mixing: Pulay every 'pulay_mixing_frequency' iterations, linear otherwise.
        
        Parameters
        ----------
        runs : int
            Current iteration number (before increment in reference code)
        preconditioner : Optional[np.ndarray]
            Preconditioner matrix, shape (n_points, n_points), default is None
            Required when use_preconditioner is True.
        
        Returns
        -------
        rho_next : np.ndarray
            Mixed density
        
        Notes
        -----
        Reference code increments runs before mixing, so we add 1 to match:
        - runs=0 (here) -> runs=1 (ref): linear (1%2=1)
        - runs=1 (here) -> runs=2 (ref): Pulay  (2%2=0)
        - runs=2 (here) -> runs=3 (ref): linear (3%2=1)
        - runs=3 (here) -> runs=4 (ref): Pulay  (4%2=0)
        """
        use_pulay = ((runs + 1) % self.pulay_mixing_frequency == 0) and (runs > 0)
        
        if runs < self.pulay_mixing_history:
            # Early phase: direct indexing
            if use_pulay:
                rho_next = self._pulay_mix_early(runs, preconditioner)
            else:
                rho_next = self._linear_mix_early(runs, preconditioner)
        else:
            # Stable phase: use rolling window (last column)
            if use_pulay:
                rho_next = self._pulay_mix_stable(preconditioner)
            else:
                rho_next = self._linear_mix_stable(preconditioner)
        
        # Correct non-positive densities (optional)
        if self.use_rho_clamp:
            negative_mask = rho_next <= 0
            rho_next[negative_mask] = self.rho_clamp_minimum
        
        return rho_next
    
    
    def _linear_mix_early(self, runs: int, preconditioner: Optional[np.ndarray]) -> np.ndarray:
        """
        Linear mixing for early phase (runs < pulay_mixing_history),
            use alternating coefficients for odd/even steps.
        
        Already incremented 'runs' before mixing, so:
        - runs=0 (here) -> runs=1 (ref, odd) -> use linear_mixing_alpha1
        - runs=1 (here) -> runs=2 (ref, even) -> use linear_mixing_alpha2
        
        Parameters
        ----------
        runs : int
            Current iteration number (before increment in reference code)
        preconditioner : Optional[np.ndarray]
            Preconditioner matrix, shape (n_points, n_points), default is None
        """

        if (runs + 1) % 2 != 0:  # equivalent to: runs % 2 == 0
            alpha = self.linear_mixing_alpha1
        else:
            alpha = self.linear_mixing_alpha2
        
        rho_in_prev  = self.rho_in_store[:, runs]
        rho_out_prev = self.rho_out_store[:, runs]
        rho_residual_prev = rho_out_prev - rho_in_prev
        
        if self.use_preconditioner:
            assert preconditioner is not None, \
                PRECONDITIONER_IS_NONE_WHEN_USE_PRECONDITIONER_IS_TRUE_ERROR
            rho_residual_prev = np.linalg.solve(preconditioner, rho_residual_prev)
        
        return rho_in_prev + alpha * rho_residual_prev
    
    
    def _linear_mix_stable(self, preconditioner: Optional[np.ndarray]) -> np.ndarray:
        """
        Linear mixing for stable phase (runs >= pulay_mixing_history),
            use alternating coefficients for odd/even steps.
        
        In the stable phase, history is managed using a rolling window.
        The previous densities are stored in the second-to-last column (-2)
        of the history buffers.
        
        Parameters
        ----------
        preconditioner : Optional[np.ndarray]
            Preconditioner matrix, shape (n_points, n_points), default is None
            Used for preconditioning if use_preconditioner is True.
        
        Returns
        -------
        rho_next : np.ndarray
            Mixed density for next iteration
        
        Notes
        -----
        Uses alternating coefficients based on iteration count:
        - Odd iterations (runs+1 % 2 != 0): use linear_mixing_alpha1
        - Even iterations (runs+1 % 2 == 0): use linear_mixing_alpha2
        """

        runs = self.iteration_count
        
        # Alternating coefficients (matching reference after increment)
        if (runs + 1) % 2 != 0:  # equivalent to: runs % 2 == 0
            alpha = self.linear_mixing_alpha1
        else:
            alpha = self.linear_mixing_alpha2
        
        rho_in_prev  = self.rho_in_store[:, -1]
        rho_out_prev = self.rho_out_store[:, -1]
        rho_residual_prev = rho_out_prev - rho_in_prev
        
        if self.use_preconditioner:
            assert preconditioner is not None, \
                PRECONDITIONER_IS_NONE_WHEN_USE_PRECONDITIONER_IS_TRUE_ERROR
            rho_residual_prev = np.linalg.solve(preconditioner, rho_residual_prev)
        
        return rho_in_prev + alpha * rho_residual_prev
    
    
    def _pulay_mix_early(self, runs: int, preconditioner: Optional[np.ndarray]) -> np.ndarray:
        """
        Pulay mixing for early phase (runs < history)

        Theoretical formula:
        --------------------
        Pulay mixing's projection matrix takes the form

            P = alpha * I - (Δρ_in + alpha*ΔF) @ (ΔF^T @ ΔF)^(-1) @ ΔF^T
              = alpha * I - (rho_in_residual + alpha * F) @ F_pinv

        where F_pinv := (F^T @ F)^(-1) @ F^T is the pseudo-inverse of F.
        
        Then, the updating formula takes the form:

            rho_out -> rho_in_prev + P @ (rho_out_prev - rho_in_prev)
                 = rho_in_prev + P @ Fk
                 = rho_in_prev + alpha * Fk - (rho_in_residual + alpha * F) @ F_pinv @ Fk
                 = rho_in_prev + alpha * Fk - alpha * F @ F_pinv @ Fk - rho_in_residual @ F_pinv @ Fk
                 = rho_in_prev + alpha * (Fk - F @ F_pinv @ Fk) - rho_in_residual @ F_pinv @ Fk
                := rho_in_prev + alpha * Fk_orthogonal - rho_in_residual @ F_pinv @ Fk
        
        where Fk_orthogonal := Fk - F @ F_pinv @ Fk is the orthogonal component of Fk to the space spanned by F.
        """

        # Compute residuals using differences
        rho_in_residual  = np.diff(self.rho_in_store[:, :runs+1], axis=1)
        rho_out_residual = np.diff(self.rho_out_store[:, :runs+1], axis=1)

        # Residual history matrix F
        F = rho_out_residual - rho_in_residual
        
        # Current residual
        rho_in_prev  = self.rho_in_store[:, runs]
        rho_out_prev = self.rho_out_store[:, runs]
        Fk = rho_out_prev - rho_in_prev

        # Pulay projection matrix formula
        F_gamma = np.linalg.inv(F.T @ F) @ (F.T @ Fk)
        Fk_orthogonal = Fk - F @ F_gamma

        if not self.use_preconditioner:
            return rho_in_prev + self.pulay_mixing_parameter * Fk_orthogonal - rho_in_residual @ F_gamma
        else:
            return rho_in_prev + self.pulay_mixing_parameter * np.linalg.solve(preconditioner, Fk_orthogonal) - rho_in_residual @ F_gamma
    
    
    def _pulay_mix_stable(self, preconditioner: Optional[np.ndarray]) -> np.ndarray:
        """
        Pulay mixing for stable phase (runs >= pulay_mixing_history).

        In the stable phase, history is managed using a rolling window.
        All history columns are used to compute the residual matrix F.

        Theoretical formula:
        --------------------
        Pulay mixing's projection matrix takes the form

            P = alpha * I - (Δρ_in + alpha*ΔF) @ (ΔF^T @ ΔF)^(-1) @ ΔF^T
              = alpha * I - (rho_in_residual + alpha * F) @ F_pinv

        where F_pinv := (F^T @ F)^(-1) @ F^T is the pseudo-inverse of F.
        
        Then, the updating formula takes the form:

            rho_out -> rho_in_prev + P @ (rho_out_prev - rho_in_prev)
                 = rho_in_prev + P @ Fk
                 = rho_in_prev + alpha * Fk - (rho_in_residual + alpha * F) @ F_pinv @ Fk
                 = rho_in_prev + alpha * Fk - alpha * F @ F_pinv @ Fk - rho_in_residual @ F_pinv @ Fk
                 = rho_in_prev + alpha * (Fk - F @ F_pinv @ Fk) - rho_in_residual @ F_pinv @ Fk
                := rho_in_prev + alpha * Fk_orthogonal - rho_in_residual @ F_pinv @ Fk
        
        where Fk_orthogonal := Fk - F @ F_pinv @ Fk is the orthogonal component of Fk to the space spanned by F.
        
        Parameters
        ----------
        preconditioner : Optional[np.ndarray]
            Preconditioner matrix, shape (n_points, n_points), default is None
            Used for preconditioning if use_preconditioner is True.
        
        Returns
        -------
        rho_next : np.ndarray
            Mixed density for next iteration
        """
        # Compute residuals using differences (all history columns)
        rho_in_residual = np.diff(self.rho_in_store, axis=1)
        rho_out_residual = np.diff(self.rho_out_store, axis=1)
        
        # Residual history matrix F
        F = rho_out_residual - rho_in_residual
        
        # Current residual
        rho_in_prev = self.rho_in_store[:, -1]
        rho_out_prev = self.rho_out_store[:, -1]
        Fk = rho_out_prev - rho_in_prev

        # Pulay projection matrix formula
        F_gamma = np.linalg.inv(F.T @ F) @ (F.T @ Fk)
        Fk_orthogonal = Fk - F @ F_gamma

        if not self.use_preconditioner:
            return rho_in_prev + self.pulay_mixing_parameter * Fk_orthogonal - rho_in_residual @ F_gamma
        else:
            return rho_in_prev + self.pulay_mixing_parameter * np.linalg.solve(preconditioner, Fk_orthogonal) - rho_in_residual @ F_gamma
    
        
