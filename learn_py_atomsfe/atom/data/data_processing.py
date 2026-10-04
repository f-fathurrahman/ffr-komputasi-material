
__author__ = "Qihao Cheng"



"""Data processing utilities for atomic DFT calculations."""

import numpy as np
from typing import Optional, Dict, Any, Literal

# Valid smooth methods
VALID_SMOOTH_METHODS = ['lowpass', 'savgol', 'moving_avg', 'spline', 'gaussian', 'exp_weighted', 'cascade']
VALID_FE_BOUNDARY_FIT_METHODS = ['spline', 'cubic_spline', 'polynomial', 'pchip', 'akima']

# Error messages
RHO_NOT_NDARRAY_ERROR = \
    "parameter 'rho' must be a numpy.ndarray, get type {} instead."
RHO_NOT_2D_ARRAY_ERROR = \
    "parameter 'rho' must be a 2D numpy.ndarray, get dimension {} instead."
VXC_TARGET_PHYSICAL_NOT_NDARRAY_ERROR = \
    "parameter 'vxc_target_physical' must be a numpy.ndarray, get type {} instead."
VXC_TARGET_PHYSICAL_NOT_2D_ARRAY_ERROR = \
    "parameter 'vxc_target_physical' must be a 2D numpy.ndarray, get dimension {} instead."
VXC_NOT_NDARRAY_ERROR = \
    "parameter 'vxc' must be a numpy.ndarray, get type {} instead."
VALUES_NOT_NDARRAY_ERROR = \
    "parameter 'values' must be a numpy.ndarray, get type {} instead."
R_NOT_NDARRAY_ERROR = \
    "parameter 'r' must be a numpy.ndarray, get type {} instead."
VXC_AND_R_NOT_SAME_LENGTH_ERROR = \
    "vxc and r must have the same length: {} vs {}."
VALUES_AND_R_NOT_SAME_LENGTH_ERROR = \
    "values and r must have the same length: {} vs {}."
R_THRESHOLD_NOT_FLOAT_ERROR = \
    "parameter 'r_threshold' must be a float, get type {} instead."
R_THRESHOLD_NOT_POSITIVE_ERROR = \
    "parameter 'r_threshold' must be positive, get {} instead."
METHOD_NOT_STRING_ERROR = \
    "parameter 'method' must be a string, get type {} instead."
METHOD_NOT_IN_VALID_LIST_ERROR = \
    "parameter 'method' must be in {}, get {} instead."
KWARGS_NOT_DICT_ERROR = \
    "parameter 'kwargs' must be a dict, get type {} instead."

UNKNOWN_SMOOTH_METHOD_ERROR = \
    "Unknown smoothing method: {}. Choose from: {}."

QUADRATURE_NODES_NOT_NDARRAY_ERROR = \
    "parameter 'quadrature_nodes' must be a numpy.ndarray, get type {} instead."
VALUES_NOT_1D_ARRAY_ERROR = \
    "parameter 'values' must be a 1D numpy.ndarray, get dimension {} instead."
QUADRATURE_NODES_NOT_1D_ARRAY_ERROR = \
    "parameter 'quadrature_nodes' must be a 1D numpy.ndarray, get dimension {} instead."
VALUES_SIZE_MISMATCH_ERROR = \
    "parameter 'values' size {} must equal n_finite_elements * n_quad_per_element = {}."
QUADRATURE_NODES_SIZE_MISMATCH_ERROR = \
    "parameter 'quadrature_nodes' size {} must equal n_finite_elements * n_quad_per_element = {}."
N_FINITE_ELEMENTS_NOT_INT_ERROR = \
    "parameter 'n_finite_elements' must be an int, get type {} instead."
N_FINITE_ELEMENTS_NOT_POSITIVE_ERROR = \
    "parameter 'n_finite_elements' must be positive, get {} instead."
N_QUAD_PER_ELEMENT_NOT_INT_ERROR = \
    "parameter 'n_quad_per_element' must be an int, get type {} instead."
N_QUAD_PER_ELEMENT_NOT_POSITIVE_ERROR = \
    "parameter 'n_quad_per_element' must be positive, get {} instead."
QUADRATURE_NODES_NOT_MONOTONIC_ERROR = \
    "parameter 'quadrature_nodes' must be strictly increasing."
N_EXTEND_N_EXCLUDE_SUM_TOO_LARGE_ERROR = \
    "n_extend + n_exclude = {} must be < n_quad_per_element = {}."


class DataProcessor:

    @staticmethod
    def _f_exp_neg_inv_x(x: np.ndarray) -> np.ndarray:
        """f(x) = exp(-1/x) for x > 0, else 0 (building block for Gevrey step)."""
        x = np.asarray(x, dtype=float)
        out = np.zeros_like(x, dtype=float)
        pos = x > 0
        out[pos] = np.exp(-1.0 / x[pos])
        return out

    @staticmethod
    def smooth_step_exp_inv_bridge(r: np.ndarray, a: float, b: float) -> np.ndarray:
        """
        S(r) = 0 for r <= a, S(r) = 1 for r >= b; on (a, b) a C^infty Gevrey step.
        Require a < b (same units as r).
        """
        if not (a < b):
            raise ValueError("Require a < b (onset a, completion b, same units as r).")
        r = np.asarray(r, dtype=float)
        fa, fb = float(a), float(b)
        out = np.ones_like(r, dtype=float)
        lo = r <= fa
        hi = r >= fb
        mid = ~(lo | hi)
        out[lo] = 0.0
        out[hi] = 1.0
        if np.any(mid):
            u = (r[mid] - fa) / (fb - fa)
            eps = 1e-15
            u = np.clip(u, eps, 1.0 - eps)
            fu = DataProcessor._f_exp_neg_inv_x(u)
            fv = DataProcessor._f_exp_neg_inv_x(1.0 - u)
            out[mid] = fu / (fu + fv)
        return out

    @staticmethod
    def smooth_rho_with_loglinear_asymptote(
        r_bohr: np.ndarray,
        rho: np.ndarray,
        r_a: float,
        r_b: float,
        r_c: float,
        *,
        rho_floor: float = 1e-30,
    ) -> np.ndarray:
        """
        Smooth radial density toward a log-linear asymptote with a C^infty blend.

        All radii in Bohr. Require ``r_a < r_b < r_c``: LS fit of log rho vs r on ``[r_a, r_c]``;
        raw rho unchanged for ``r <= r_b``; blend to ``rho_asym = exp(alpha + beta r)`` on ``[r_b, r_c]``;
        pure asymptote for ``r >= r_c``.
        """
        if not (r_a < r_b < r_c):
            raise ValueError("Require r_a < r_b < r_c (Bohr).")
        r = np.asarray(r_bohr, dtype=float)
        y = np.asarray(rho, dtype=float)
        if r.shape != y.shape:
            raise ValueError("r_bohr and rho must have the same shape.")

        fit = (r >= r_a) & (r <= r_c)
        if not np.any(fit):
            raise ValueError("No grid points fall in [r_a, r_c] for the fit.")
        log_rho = np.log(np.maximum(y, float(rho_floor)))
        A = np.column_stack([np.ones(int(np.count_nonzero(fit))), r[fit]])
        bvec = log_rho[fit]
        coeff, _, rank, _ = np.linalg.lstsq(A, bvec, rcond=None)
        if rank < 2:
            raise ValueError(
                "Least-squares fit is rank-deficient; need at least two distinct r in the fit window."
            )

        alpha, beta = float(coeff[0]), float(coeff[1])
        rho_asym = np.exp(alpha + beta * r)

        w = DataProcessor.smooth_step_exp_inv_bridge(r, r_b, r_c)
        rho_smooth = (1.0 - w) * y + w * rho_asym
        return rho_smooth


    @staticmethod
    def smooth_radial_data(
        values      : np.ndarray,
        r           : np.ndarray,
        r_threshold : float = 10.0,
        method      : str   = 'savgol',
        **kwargs
    ):
        """
        Smooth radial values (e.g. vxc, exc) for r > r_threshold to reduce numerical instability and high-frequency oscillations.

        Common smoothing methods for filtering high-frequency noise:
        1. 'lowpass'      - Low-pass Butterworth filter (RECOMMENDED for high-frequency noise)
        2. 'savgol'       - Savitzky-Golay filter : a local polynomial regression-based smoothing method that reduces
                            noise while preserving the shape, height, and position of signal features such as peaks,
                            making it particularly suitable for spectroscopic and physical data analysis.

        3. 'moving_avg'   - Moving average (simple, good for high-frequency filtering)
        4. 'spline'       - Spline smoothing (controllable smoothness)
        5. 'gaussian'     - Gaussian filter (good smoothing, adjustable strength)
        6. 'exp_weighted' - Exponentially weighted moving average (smooth for large r)
        7. 'cascade'      - Apply multiple smoothing methods in sequence (strongest filtering)


        Parameters
        ----------
        values : np.ndarray
            Radial values to smooth (e.g. vxc, exc)
        r : np.ndarray
            Radius values (same length as values)
        r_threshold : float
            Threshold radius. Values with r > r_threshold will be smoothed.
        method : str
            Smoothing method: 'lowpass' , 'savgol'(default), 'moving_avg', 'spline', 'gaussian', 'exp_weighted', 'cascade'
        **kwargs : Optional[Dict]
            Additional parameters for smoothing methods:
            - lowpass: cutoff (default: 0.05), order (default: 6) - lower cutoff = stronger filtering
            - savgol: window_length (default: min(30% of data, len(data)//2*2+1)), polyorder (default: 2)
            - moving_avg: window_size (default: 25% of data, min 25) - larger = stronger filtering
            - spline: s (smoothing factor, default: len(data) * variance * 0.8) - larger = stronger
            - gaussian: sigma (default: max(2.0, 1% of data length)) - larger = stronger filtering
            - exp_weighted: alpha (default: 0.15) - smaller = stronger filtering
            - cascade: methods (list), kwargs_list (list of kwargs for each method)

        Returns
        -------
        values_smoothed : np.ndarray
            Smoothed values (only r > r_threshold is smoothed)
        """

        # Convert to numpy arrays
        values = np.asarray(values).copy()
        r = np.asarray(r)

        # Check arguments
        assert isinstance(values, np.ndarray), \
            VALUES_NOT_NDARRAY_ERROR.format(type(values))
        assert isinstance(r, np.ndarray), \
            R_NOT_NDARRAY_ERROR.format(type(r))
        assert len(values) == len(r), \
            VALUES_AND_R_NOT_SAME_LENGTH_ERROR.format(len(values), len(r))
        assert isinstance(r_threshold, float), \
            R_THRESHOLD_NOT_FLOAT_ERROR.format(type(r_threshold))
        assert r_threshold > 0, \
            R_THRESHOLD_NOT_POSITIVE_ERROR.format(r_threshold)
        assert isinstance(method, str), \
            METHOD_NOT_STRING_ERROR.format(type(method))
        assert method in VALID_SMOOTH_METHODS, \
            METHOD_NOT_IN_VALID_LIST_ERROR.format(VALID_SMOOTH_METHODS, method)
        assert isinstance(kwargs, dict), \
            KWARGS_NOT_DICT_ERROR.format(type(kwargs))
        
        # Find indices where r > r_threshold
        large_r_mask = r > r_threshold
        
        if not np.any(large_r_mask):
            # No points to smooth
            return values

        # Extract data for smoothing
        values_large_r = values[large_r_mask]
        r_large_r      = r[large_r_mask]

        if len(values_large_r) < 3:
            # Not enough points to smooth
            return values
        
        # Apply smoothing based on method
        if method == 'lowpass':
            # Low-pass Butterworth filter - BEST for filtering high-frequency oscillations
            from scipy.signal import butter, filtfilt
            # Normalized cutoff frequency (0 < cutoff < 1, where 1 is Nyquist frequency)
            # Lower cutoff = stronger filtering (removes more high frequencies)
            cutoff = kwargs.get('cutoff', 0.05)  # Default: filter out frequencies > 5% of Nyquist (stronger)
            order = kwargs.get('order', 6)  # Filter order (higher = sharper cutoff, increased default)
            
            # Design Butterworth low-pass filter
            b, a = butter(order, cutoff, btype='low', analog=False)
            
            # Apply filter (filtfilt applies forward and backward for zero phase distortion)
            values_smoothed_large_r = filtfilt(b, a, values_large_r)

        elif method == 'savgol':
            from scipy.signal import savgol_filter
            # Use larger window and lower polyorder for better high-frequency filtering
            # Default: use up to 30% of data points for window (much larger for stronger filtering)
            default_window = min(int(len(values_large_r) * 0.3) // 2 * 2 + 1, len(values_large_r))
            default_window = max(default_window, 21)  # At least 21
            window_length = kwargs.get('window_length', default_window)
            polyorder = kwargs.get('polyorder', 2)  # Lower order for stronger smoothing
            # Ensure window_length is odd and <= len(data)
            if window_length % 2 == 0:
                window_length += 1
            window_length = min(window_length, len(values_large_r))
            if window_length < polyorder + 2:
                polyorder = max(1, window_length - 2)

            values_smoothed_large_r = savgol_filter(values_large_r, window_length, polyorder)
        
        elif method == 'moving_avg':
            # Much larger window for better high-frequency filtering
            # Default: use up to 25% of data points (much stronger filtering)
            default_window = max(int(len(values_large_r) * 0.25), 25)  # At least 25, up to 25% of data
            window_size = kwargs.get('window_size', default_window)
            window_size = min(window_size, len(values_large_r))
            # Use convolution for moving average
            kernel = np.ones(window_size) / window_size
            values_smoothed_large_r = np.convolve(values_large_r, kernel, mode='same')
            # Handle boundaries
            for i in range(window_size // 2):
                values_smoothed_large_r[i] = np.mean(values_large_r[:i+window_size//2+1])
                values_smoothed_large_r[-(i+1)] = np.mean(values_large_r[-(i+window_size//2+1):])
        
        elif method == 'spline':
            from scipy.interpolate import UnivariateSpline
            # Calculate smoothing factor based on data variance
            # Larger s = stronger smoothing (filters more high-frequency noise)
            data_variance = np.var(values_large_r)
            # Default: use 80% of variance (much stronger smoothing)
            s = kwargs.get('s', len(values_large_r) * data_variance * 0.8)
            spline = UnivariateSpline(r_large_r, values_large_r, s=s, k=3)
            values_smoothed_large_r = spline(r_large_r)
        
        elif method == 'gaussian':
            from scipy.ndimage import gaussian_filter1d
            # Larger sigma = stronger smoothing (better for high-frequency filtering)
            # Default: adaptive sigma based on data length (stronger for longer sequences)
            default_sigma = max(2.0, len(values_large_r) * 0.01)  # At least 2.0, or 1% of data length
            sigma = kwargs.get('sigma', default_sigma)
            values_smoothed_large_r = gaussian_filter1d(values_large_r, sigma=sigma)
        
        elif method == 'exp_weighted':
            # Lower alpha = stronger smoothing (more weight on previous values)
            alpha = kwargs.get('alpha', 0.15)  # Further decreased for stronger filtering
            values_smoothed_large_r = np.zeros_like(values_large_r)
            values_smoothed_large_r[0] = values_large_r[0]
            for i in range(1, len(values_large_r)):
                values_smoothed_large_r[i] = alpha * values_large_r[i] + (1 - alpha) * values_smoothed_large_r[i-1]
        
        elif method == 'cascade':
            # Apply multiple smoothing methods in sequence for strongest filtering
            methods_list = kwargs.get('methods', ['lowpass', 'moving_avg'])
            # Default: stronger parameters for cascade
            default_window = max(int(len(values_large_r) * 0.2), 20)
            kwargs_list = kwargs.get('kwargs_list', [
                {'cutoff': 0.05, 'order': 6},  # Strong lowpass
                {'window_size': default_window}  # Large moving average
            ])

            values_smoothed_large_r = values_large_r.copy()
            for m, kws in zip(methods_list, kwargs_list):
                # Apply each smoothing method directly (avoid recursion)
                if m == 'lowpass':
                    from scipy.signal import butter, filtfilt
                    cutoff = kws.get('cutoff', 0.05)
                    order = kws.get('order', 6)
                    b, a = butter(order, cutoff, btype='low', analog=False)
                    values_smoothed_large_r = filtfilt(b, a, values_smoothed_large_r)
                elif m == 'moving_avg':
                    window_size = kws.get('window_size', default_window)
                    window_size = min(window_size, len(values_smoothed_large_r))
                    kernel = np.ones(window_size) / window_size
                    values_smoothed_large_r = np.convolve(values_smoothed_large_r, kernel, mode='same')
                elif m == 'gaussian':
                    from scipy.ndimage import gaussian_filter1d
                    sigma = kws.get('sigma', max(2.0, len(values_smoothed_large_r) * 0.01))
                    values_smoothed_large_r = gaussian_filter1d(values_smoothed_large_r, sigma=sigma)
                else:
                    # For other methods, use recursive call (but with r_threshold=0 to smooth all)
                    values_smoothed_large_r = DataProcessor.smooth_radial_data(
                        values_smoothed_large_r, r_large_r,
                        r_threshold=-1.0,  # Negative threshold to smooth all points
                        method=m, **kws
                    )
        
        else:
            raise ValueError(UNKNOWN_SMOOTH_METHOD_ERROR.format(method, VALID_SMOOTH_METHODS))
        
        # Replace smoothed values
        values[large_r_mask] = values_smoothed_large_r

        return values

    # Backward compatibility alias
    smooth_vxc_data = smooth_radial_data


    @staticmethod
    def smooth_radial_data_at_finite_elements_boundaries(
        values             : np.ndarray,
        quadrature_nodes   : np.ndarray,
        n_finite_elements  : int,
        n_quad_per_element : int,
        n_extend           : int = 7,
        n_exclude          : int = 3,
        r_threshold        : float = 10.0,
    ) -> np.ndarray:
        """
        Smooth radial values at finite element boundaries to reduce discontinuities.

        Assumes data follows exponential growth/decay: |v| ~ exp(a*r + b).
        Fits log(|v| + eps) = a*r + b in least-squares sense, then replaces the
        target interval with sign(v_original) * exp(a*r + b).

        For each boundary (index i = last of elem k, i+1 = first of elem k+1):
        - Fit using left points [i-n_extend, i-n_exclude] and right [i+1+n_exclude, i+1+n_extend]
        - Replace values at [i-n_exclude+1, i+1+n_exclude-1] with the fitted curve

        Parameters
        ----------
        values : np.ndarray
            1D radial values to smooth.
        quadrature_nodes : np.ndarray
            1D quadrature node positions (must be strictly increasing).
        n_finite_elements : int
            Number of finite elements.
        n_quad_per_element : int
            Quadrature points per element.
        n_extend : int, default=7
            Number of points to use for fitting on each side of the boundary.
        n_exclude : int, default=3
            Number of points to exclude from fitting (closest to boundary).
        r_threshold : float, default=10.0
            Only smooth boundaries where r > r_threshold.

        Returns
        -------
        np.ndarray
            Smoothed values (copy).
        """
        raise DeprecationWarning("Function 'smooth_radial_data_at_finite_elements_boundaries' is now deprecated. This method is not helpful for radial data smoothing.")

        # Type check
        assert isinstance(values, np.ndarray), \
            VALUES_NOT_NDARRAY_ERROR.format(type(values))
        assert values.ndim == 1, \
            VALUES_NOT_1D_ARRAY_ERROR.format(values.ndim)
        assert isinstance(quadrature_nodes, np.ndarray), \
            QUADRATURE_NODES_NOT_NDARRAY_ERROR.format(type(quadrature_nodes))
        assert quadrature_nodes.ndim == 1, \
            QUADRATURE_NODES_NOT_1D_ARRAY_ERROR.format(quadrature_nodes.ndim)
        assert isinstance(r_threshold, float), \
            R_THRESHOLD_NOT_FLOAT_ERROR.format(type(r_threshold))
        assert r_threshold > 0, \
            R_THRESHOLD_NOT_POSITIVE_ERROR.format(r_threshold)
        assert isinstance(n_finite_elements, int), \
            N_FINITE_ELEMENTS_NOT_INT_ERROR.format(type(n_finite_elements))
        assert n_finite_elements > 0, \
            N_FINITE_ELEMENTS_NOT_POSITIVE_ERROR.format(n_finite_elements)
        assert isinstance(n_quad_per_element, int), \
            N_QUAD_PER_ELEMENT_NOT_INT_ERROR.format(type(n_quad_per_element))
        assert n_quad_per_element > 0, \
            N_QUAD_PER_ELEMENT_NOT_POSITIVE_ERROR.format(n_quad_per_element)
        assert isinstance(n_extend, int) and n_extend > 0, \
            "n_extend must be a positive int."
        assert isinstance(n_exclude, int) and n_exclude > 0, \
            "n_exclude must be a positive int."
        assert n_extend + n_exclude < n_quad_per_element, \
            N_EXTEND_N_EXCLUDE_SUM_TOO_LARGE_ERROR.format(n_extend + n_exclude, n_quad_per_element)

        # Size check
        expected_size = n_finite_elements * n_quad_per_element
        assert values.size == expected_size, \
            VALUES_SIZE_MISMATCH_ERROR.format(values.size, expected_size)
        assert quadrature_nodes.size == expected_size, \
            QUADRATURE_NODES_SIZE_MISMATCH_ERROR.format(quadrature_nodes.size, expected_size)

        # Monotonicity check
        assert np.all(np.diff(quadrature_nodes) > 0), QUADRATURE_NODES_NOT_MONOTONIC_ERROR

        # Reshape: boundaries are between elements (last of row k, first of row k+1)
        v_out = values.copy()

        # Boundary indices: i = last index of element k, i+1 = first of element k+1
        for k in range(n_finite_elements - 1):
            i = (k + 1) * n_quad_per_element - 1
            r_boundary = 0.5 * (quadrature_nodes[i] + quadrature_nodes[i + 1])
            if r_boundary <= r_threshold:
                continue

            # Fit range: left [i-n_extend, i-n_exclude], right [i+1+n_exclude, i+1+n_extend]
            left_fit_idx = np.arange(i - n_extend, i - n_exclude + 1)
            right_fit_idx = np.arange(i + 1 + n_exclude, i + 1 + n_extend + 1)
            fit_idx = np.concatenate([left_fit_idx, right_fit_idx])
            r_fit = quadrature_nodes[fit_idx]
            v_fit = v_out[fit_idx]

            # Replace range: [i-n_exclude+1, i+1+n_exclude-1] = [i-2, i-1, i, i+1, i+2, i+3] for n_exclude=3
            replace_idx = np.arange(i - n_exclude + 1, i + n_exclude + 1)
            r_replace = quadrature_nodes[replace_idx]

            # Exponential fit in log space: log(|v| + eps) = a*r + b
            eps = 1e-30
            log_v_fit = np.log(np.abs(v_fit) + eps)
            coeffs = np.polyfit(r_fit, log_v_fit, deg=1)
            a, b = coeffs[0], coeffs[1]
            v_replaced = np.exp(a * r_replace + b)
            v_out[replace_idx] = np.sign(v_out[replace_idx]) * v_replaced

        return v_out




    @staticmethod
    def symlog(x, linthresh=0.002):
        """
        Symmetric logarithm transformation.
        
        - |x| <= linthresh: y = x / linthresh (linear)
        - |x| > linthresh:  y = sign(x) * [log(|x|/linthresh) + 1] (logarithmic)
        
        Parameters
        ----------
        x : array-like
            Input values
        linthresh : float
            Linear threshold for symlog
        
        Returns
        -------
        Transformed values
        """
        x = np.asarray(x, dtype=float)
        abs_x = np.abs(x)
        out = np.empty_like(x)
        linear = abs_x <= linthresh
        out[linear] = x[linear] / linthresh
        out[~linear] = np.sign(x[~linear]) * (np.log(abs_x[~linear] / linthresh) + 1)
        return out


    @staticmethod
    def symexp(y, linthresh=0.002):
        """
        Inverse of symlog (symmetric exponential).
        
        - |y| <= 1: x = y * linthresh (linear)
        - |y| > 1:  x = sign(y) * linthresh * exp(|y| - 1) (exponential)
        
        Parameters
        ----------
        y : array-like
            Symlog-transformed values
        linthresh : float
            Linear threshold (must match symlog)
        
        Returns
        -------
        Original values
        """
        y = np.asarray(y, dtype=float)
        abs_y = np.abs(y)
        out = np.empty_like(y)
        linear = abs_y <= 1
        out[linear] = y[linear] * linthresh
        out[~linear] = np.sign(y[~linear]) * linthresh * np.exp(abs_y[~linear] - 1)
        return out


    @staticmethod
    def symlog_derivative(x, linthresh=0.002):
        """
        Absolute value of symlog derivative: |d/dx symlog(x)|.
        
        Used for weight calculation when targets are symlog-transformed.
        - For |x| <= linthresh (linear region): |symlog'(x)| = 1 / linthresh
        - For |x| > linthresh (log region): |symlog'(x)| = 1 / |x|
        
        Parameters
        ----------
        x : array-like
            Input values (in physical space, before symlog)
        linthresh : float
            Linear threshold (must match symlog)
        
        Returns
        -------
        np.ndarray
            |symlog'(x)|, same shape as x
        """
        x = np.asarray(x)
        abs_x = np.abs(x)
        result = np.empty_like(abs_x, dtype=float)
        small_mask = abs_x <= linthresh
        result[small_mask] = 1.0 / linthresh
        result[~small_mask] = 1.0 / abs_x[~small_mask]
        return result


    # Function to calculate correct weights for symlog-transformed targets
    @staticmethod
    def calculate_symlog_weights(rho, vxc_target_physical, linthresh, min_weight=1e-6):
        """
        Calculate weights for loss function when targets are symlog-transformed.
        
        The weight is designed to approximate the original rho-weighted loss in physical space:
        L_original = rho * |vxc_pred - vxc_true|
        
        IMPORTANT: This function uses the target Vxc (not delta_vxc) for weight calculation,
        as the weight should be based on the target potential magnitude.
        
        Parameters
        ----------
        rho : array-like
            Density values (physical space)
        vxc_target_physical : array-like
            Target vxc values in physical space (v_x + v_c)
            This is used instead of delta_vxc because weights should reflect
            the magnitude of the target potential
        linthresh : float
            Linear threshold for symlog (default: 0.001)
        
        Returns
        -------
        weights : array
            Weights for loss function
        """
        rho = np.asarray(rho)
        vxc_target_physical = np.asarray(vxc_target_physical)

        # Type check
        assert isinstance(rho, np.ndarray), \
            RHO_NOT_NDARRAY_ERROR.format(type(rho))
        assert isinstance(vxc_target_physical, np.ndarray), \
            VXC_TARGET_PHYSICAL_NOT_NDARRAY_ERROR.format(type(vxc_target_physical))

        # Shape check
        if rho.ndim == 1:
            rho = rho.reshape(-1, 1)
        if vxc_target_physical.ndim == 1:
            vxc_target_physical = vxc_target_physical.reshape(-1, 1)
        
        # Dimension check
        assert rho.ndim == 2, \
            RHO_NOT_2D_ARRAY_ERROR.format(type(rho))
        assert vxc_target_physical.ndim == 2, \
            VXC_TARGET_PHYSICAL_NOT_2D_ARRAY_ERROR.format(type(vxc_target_physical))


        # Calculate absolute value of symlog derivative using target vxc
        abs_symlog_derivative = DataProcessor.symlog_derivative(
            vxc_target_physical, linthresh=linthresh
        )
        
        # Step 3: Calculate weight = rho / |symlog'(vxc)|
        # For |vxc| <= linthresh: weight = rho * linthresh
        # For |vxc| > linthresh: weight = rho * |vxc|
        weights = rho / abs_symlog_derivative
        
        # Step 4: Ensure weights are non-negative
        # Take absolute value of rho to ensure non-negative weights
        # (rho should already be non-negative, but this is a safety check)
        weights = np.abs(weights)
        
        # Additional safety: clip any negative or zero weights to a small positive value
        # This should not happen if rho is positive, but just in case
        weights = np.maximum(weights, min_weight)
        
        return weights



    @staticmethod
    def normalize_weights_by_atom(weights, atomic_numbers, min_weight_ratio=1e-2):
        """
        Normalize weights so that each atom contributes equally to the total loss.
        
        This is important because different atoms have different numbers of electrons
        and different numbers of data points. Without normalization, atoms with more
        data points or higher electron density would dominate the loss.
        
        The normalization ensures:
        - Each atom's total weighted contribution is equal
        - Atoms with more data points will have lower per-sample weights
        - Atoms with fewer data points will have higher per-sample weights
        - All weights are bounded below by max(weights) * min_weight_ratio to avoid numerical issues
        
        Parameters
        ----------
        weights : array-like
            Original weights for each sample
        atomic_numbers : array-like
            Atomic numbers for each sample (same length as weights)
        min_weight_ratio : float
            Minimum weight ratio relative to maximum weight (default: 1e-2).
            The minimum weight will be max(normalized_weights) * min_weight_ratio.
            All normalized weights will be at least this value to avoid numerical issues.
        
        Returns
        -------
        normalized_weights : array
            Normalized weights such that each atom contributes equally,
            with all weights >= max(normalized_weights) * min_weight_ratio
        """
        weights = np.asarray(weights)
        atomic_numbers = np.asarray(atomic_numbers)
        
        if len(weights) != len(atomic_numbers):
            raise ValueError(f"weights and atomic_numbers must have the same length: "
                            f"{len(weights)} vs {len(atomic_numbers)}")
        
        # Get unique atoms
        unique_atoms = np.unique(atomic_numbers)
        
        # Calculate total weight per atom
        atom_total_weights = {}
        atom_counts = {}
        for atom_z in unique_atoms:
            atom_mask = atomic_numbers == atom_z
            atom_total_weights[atom_z] = np.sum(weights[atom_mask])
            atom_counts[atom_z] = np.sum(atom_mask)
        
        # Calculate the target total weight per atom (use mean of all atoms)
        target_total_weight = np.mean(list(atom_total_weights.values()))
        
        # Normalize weights for each atom
        normalized_weights = weights.copy()
        for atom_z in unique_atoms:
            atom_mask = atomic_numbers == atom_z
            atom_total = atom_total_weights[atom_z]
            
            if atom_total > 0:
                # Scale weights so that this atom's total equals the target
                # normalized_weight = original_weight * (target_total / atom_total)
                normalized_weights[atom_mask] = weights[atom_mask] * (target_total_weight / atom_total)
            else:
                # If atom has zero total weight, set to uniform (shouldn't happen in practice)
                # Use a small positive value to avoid division issues
                normalized_weights[atom_mask] = target_total_weight / max(atom_counts[atom_z], 1)
        
        # Apply minimum weight bound after normalization
        # Minimum weight is max(normalized_weights) * min_weight_ratio
        max_weight = np.max(normalized_weights)
        min_weight = max_weight * min_weight_ratio
        normalized_weights = np.maximum(normalized_weights, min_weight)
        
        return normalized_weights



    @classmethod
    def inverse_transform_features(cls, X_transformed, scaler_X=None, transform_params=None, feature_idx=None):
        """
        Inverse transform features back to physical values.
        
        Parameters
        ----------
        X_transformed : array-like
            Transformed features (after scaling and/or symlog)
            Can be a single feature column or full feature matrix
        scaler_X : RobustScaler or None
            Scaler used for features (if any)
        transform_params : dict or None
            Transform parameters from prepare_data
        feature_idx : int or None
            If provided and X_transformed is a single column, specifies which feature index
            this column corresponds to (0-based). If None, assumes it's the first feature.
            
        Returns
        -------
        X_physical : array
            Features in physical space (same shape as input)
        """
        X = np.asarray(X_transformed).copy()
        original_shape = X.shape
        
        # Handle case where X is a single column but scaler expects full feature matrix
        if scaler_X is not None:
            # Check if X has fewer columns than scaler expects
            if X.ndim == 1:
                X = X.reshape(-1, 1)
            n_features_input = X.shape[1]
            n_features_scaler = scaler_X.n_features_in_ if hasattr(scaler_X, 'n_features_in_') else None
            
            if n_features_scaler is not None and n_features_input < n_features_scaler:
                # Need to create a full feature matrix for inverse transform
                # Create a matrix with zeros for other features
                X_full = np.zeros((X.shape[0], n_features_scaler))
                if feature_idx is None:
                    feature_idx = 0
                X_full[:, feature_idx] = X[:, 0]
                X = X_full
            elif n_features_input == 1 and feature_idx is not None:
                # Single column with known index - create full matrix
                if n_features_scaler is not None:
                    X_full = np.zeros((X.shape[0], n_features_scaler))
                    X_full[:, feature_idx] = X[:, 0]
                    X = X_full
        
        # Step 1: Inverse scaling (if scaler was used)
        if scaler_X is not None:
            X = scaler_X.inverse_transform(X)
        
        # Step 2: Inverse symlog (if symlog was used)
        if transform_params is not None and transform_params.get('use_symlog_features', False):
            linthresh = transform_params.get('linthresh', 0.002)
            X = cls.symexp(X, linthresh=linthresh)
        
        # If we expanded to full matrix, extract only the requested feature
        if feature_idx is not None and X.shape[1] > 1:
            X = X[:, feature_idx:feature_idx+1]
        
        # Restore original shape if input was 1D
        if len(original_shape) == 1 and X.ndim == 2 and X.shape[1] == 1:
            X = X.flatten()
        
        return X



    @classmethod
    def inverse_transform_predictions(
        cls,
        y_pred           : np.ndarray,
        scaler_y         : Optional[Any] = None,
        transform_params : Optional[Dict] = None
    ) -> np.ndarray:
        """
        Inverse transform predictions back to physical space.
        
        Parameters
        ----------
        y_pred : np.ndarray
            Predictions (transformed space)
        scaler_y : Scaler or None
            Target scaler (if used)
        transform_params : dict or None
            Transform parameters
        
        Returns
        -------
        np.ndarray
            Predictions in physical space
        """
        y = np.asarray(y_pred).copy()
        
        # Step 1: Inverse scaling
        if scaler_y is not None:
            if y.ndim == 1:
                y = y.reshape(-1, 1)
            y = scaler_y.inverse_transform(y)
        
        # Step 2: Inverse symlog
        if transform_params is not None and transform_params.get('use_symlog_targets', False):
            linthresh = transform_params.get('linthresh', 0.002)
            y = cls.symexp(y, linthresh=linthresh)
        
        return y.flatten() if y.ndim > 1 and y.shape[1] == 1 else y



