from dataclasses import dataclass

import numpy as np
from typing import Dict, Any, Optional, Tuple
import warnings

MATRIX_ASSEMBLY_DTYPE = np.float64


class Quadrature1D:
    """Factory for Gauss-Legendre / Lobatto nodes and weights on [-1, 1]."""

    @staticmethod
    def gauss_legendre(n: int) -> Tuple[np.ndarray, np.ndarray]:
        r"""
        Gauss-Legendre quadrature with Golub-Welsch algorithm
            \int_{-1}^1 f(x) dx ≈ \sum_{i=1}^n w_i f(x_i)
        Return nodes (x_i) and weights (w_i) for Gauss-Legendre on [-1,1].

        Ref: 
        - https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature
        - https://en.wikipedia.org/wiki/Gaussian_quadrature#The_Golub-Welsch_algorithm

        Note: 
        - Exact for polynomials up to degree 2n-1
        - Does not include the boundaries
        """
        assert isinstance(n, int)
        assert n > 0

        beta = 0.5 / np.sqrt(1.0 - 1.0 / (2.0 * (np.arange(1, n)))**2)
        T = np.diag(beta, k=1) + np.diag(beta, k=-1)
        nodes, eigvecs = np.linalg.eigh(T, UPLO='L')
        index = np.arange(nodes.size)
        legendre_weights = 2*(eigvecs[0, index])**2
        return nodes, legendre_weights


    
    @staticmethod
    def gauss_legendre_on_interval(
        n       : int,
        min_val : float = -1.0,
        max_val : float = 1.0,
        tol     : float = 1e-13,
    ) -> Tuple[np.ndarray, np.ndarray]:
        r"""
        Gauss-Legendre quadrature on a general interval [min_val, max_val] using
        classic Newton iteration for Legendre polynomial roots.

            \int_{min_val}^{max_val} f(x) dx ≈ \sum_{i=1}^n w_i f(x_i)

        Parameters
        ----------
        n : int
            Number of quadrature points.
        min_val, max_val : float
            Interval endpoints (min_val < max_val).
        tol : float
            Tolerance for Newton iteration on the Legendre roots.

        Returns
        -------
        mesh : np.ndarray
            Quadrature nodes on [min_val, max_val].
        weights : np.ndarray
            Corresponding quadrature weights.
        """

        # basic type / value checks (keep style consistent with other quadrature routines)
        assert isinstance(n, int)
        assert n > 0
        if not isinstance(min_val, (int, float)):
            raise TypeError("INTERVAL_ENDPOINT_NOT_REAL_NUMBER_ERROR")
        if not isinstance(max_val, (int, float)):
            raise TypeError("INTERVAL_ENDPOINT_NOT_REAL_NUMBER_ERROR")
        if not (max_val > min_val):
            raise RuntimeError("INTERVAL_ENDPOINT_MAX_NOT_GREATER_THAN_MIN_ERROR")
        if tol <= 0:
            raise RuntimeError("TOLERANCE_NOT_POSITIVE_ERROR")

        # Keep internal variable names and algorithm consistent with legacy implementation
        tolerance = tol
        number = n

        length = (max_val - min_val) / 2
        mean = (max_val + min_val) / 2
        mesh = np.zeros(number)
        weights = np.zeros(number)
        for index in range(1, int(np.floor((number + 1) / 2))+1):
            z = np.cos(np.pi * (index - 0.25) / (number + 0.5))
            flag = True
            while flag:
                p1 = 1.0
                p2 = 0.0
                for j in range(1, number + 1):
                    p3 = p2
                    p2 = p1
                    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
                pp = number * (p2 - z * p1) / (1.0 - z * z)
                z1 = z
                z = z1 - p1 / pp
                flag = (abs(z - z1) > tolerance)
            mesh[index - 1] = mean - length * z
            mesh[number - index] = mean + length * z
            weights[index - 1] = 2 * length / ((1.0 - z * z) * (pp * pp))
            weights[number - index] = weights[index - 1]
        return mesh, weights



    @staticmethod
    def lobatto(n: int, tol: float = 1e-14) -> Tuple[np.ndarray, np.ndarray]:
        r"""
        Legendre-Gauss-lobatto points with Chebyshev-Lobatto initial guess
            \int_{-1}^1 f(x) dx ≈ \sum_{i=1}^n w_i f(x_i)
        Return nodes (x_i) and weights (w_i) for Legendre-Gauss-Lobatto on [-1,1].
        
        Ref: 
        - Greg von Winckel (2024). Legende-Gauss-Lobatto nodes and weights
        - https://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights

        Note:
        - Exact for polynomials up to degree 2n-3
        - Includes the boundaries
        """

        assert isinstance(n, int)
        assert n > 0

        # initial guess: Chebyshev–Lobatto
        nodes = np.cos(np.pi * np.arange(n+1) / n)
        nodes_old = 2*np.ones_like(nodes)
        Pvals = np.zeros((n+1, n+1))
        while np.max(np.abs(nodes - nodes_old)) > tol:
            nodes_old = nodes
            Pvals[0, :] = 1.0
            Pvals[1, :] = nodes
            for k in range(2, n+1):
                Pvals[k, :] = ((2*k-1)*nodes*Pvals[k-1,:] - (k-1)*Pvals[k-2,:]) / k
            nodes = nodes_old - (nodes*Pvals[n,:] - Pvals[n-1,:]) / ((n+1)*Pvals[n,:])
        legendre_weights = 2.0 / (n*(n+1) * (Pvals[n,:]**2))
        return nodes, legendre_weights



class Mesh1D:
    """1D radial mesh on [0, domain_size] with uniform / polynomial / exponential spacing."""

    def __init__(
        self,
        domain_size            : float,                 # extent of radial domain, in Bohr
        finite_elements_number : int,                   # number of finite elements
        mesh_type              : str,                   # 'uniform', 'polynomial', 'exponential'
        clustering_param       : float | None,          # concentration parameter 'for' node distribution
        exp_shift              : float | None  = None,  # shift parameter 'for' exponential distribution

        # Deprecated parameters
        finite_elements_num    : Optional[int] = None,  # Deprecated: use finite_elements_number instead
    ):

        if finite_elements_num is not None:
            if finite_elements_number is not None:
                raise RuntimeError("Both 'finite_elements_num' and 'finite_elements_number' are provided, please only use one of them.")
            finite_elements_number = finite_elements_num
            # Decided not to raise warnings here
            # print(FINITE_ELEMENTS_NUM_DEPRECATED_WARNING.format(finite_elements_num))

        self.domain_size   = float(domain_size)
        self.n_elem        = finite_elements_number
        self.mesh_type     = mesh_type
        self.concentration = clustering_param
        self.exp_shift     = exp_shift
        self.check_initial_parameters()


    def check_initial_parameters(self):
        # domain size
        assert isinstance(self.domain_size, float)
        assert self.domain_size > 0.0
        
        # number of finite elements
        assert isinstance(self.n_elem, int)
        assert self.n_elem > 0
        
        # mesh type
        assert isinstance(self.mesh_type, str)
        assert self.mesh_type in ['uniform', 'polynomial', 'exponential']
        
        # clustering parameter 
        if self.mesh_type in ['polynomial', 'exponential']:
            assert isinstance(self.concentration, float)
            assert self.concentration > 0.
            if self.mesh_type == 'exponential' and self.concentration == 1.0:
                raise NotImplementedError("Exponential mesh with concentration 1.0 is not implemented yet, should degrade to uniform mesh")
        elif self.mesh_type == 'uniform':
            if self.concentration is not None:
                print("CLUSTERING_PARAMETER_NOT_NONE_FOR_UNIFORM_MESH_TYPE_WARNING")
        else:
            raise RuntimeError("This error should never be raised.")

        # exp_shift parameter 
        if self.mesh_type == 'exponential':
            if self.exp_shift is None:
                self.exp_shift = 0.0  # default value
            assert isinstance(self.exp_shift, float)
            assert self.exp_shift >= 0.
        elif self.mesh_type in ['uniform', 'polynomial']:
            if self.exp_shift is not None:
                print("EXP_SHIFT_NOT_NONE_FOR_UNIFORM_AND_POLYNOMIAL_MESH_TYPE_WARNING")
        else:
            raise RuntimeError("This error should never be raised.")


    def generate_mesh_nodes_and_width(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate mesh nodes in [0, domain_size] for uniform, polynomial, and exponential distributions.
        For polynomial and exponential distributions, the nodes are distributed according to the clustering parameter (s or a) and exp_shift (c).

        Mesh types:
        - Uniform
            r_i = i * domain_size / n_elem
        - Polynomial:
            r_i = domain_size * (i**s) / (n_elem ** s)
        - Exponential:
            r_i = c + (domain_size-c) * (a**(i/(n_elem-1)) - 1) / (a**(n_elem/(n_elem-1)) - 1)
        
        Returns:
            mesh_nodes: array of mesh nodes r in [0, domain_size]
            mesh_width: array of mesh widths
        """
        if self.mesh_type == "uniform":
            mesh_nodes = np.linspace(0.0, self.domain_size, self.n_elem + 1)
        elif self.mesh_type == "polynomial":
            s = float(self.concentration)
            mesh_nodes = self.domain_size * (np.arange(self.n_elem+1)**s) / (self.n_elem ** s)
        elif self.mesh_type == "exponential":
            a = float(self.concentration)
            c = float(self.exp_shift)
            beta = np.log(a) / (self.n_elem - 1)
            alpha = (self.domain_size - c) / (np.exp(beta*self.n_elem) - 1.0)
            mesh_nodes = c + alpha * (np.exp(beta*np.arange(self.n_elem + 1)) - 1.0)
            mesh_nodes[0]  = 0.0 
            mesh_nodes[-1] = self.domain_size
        else:
            raise RuntimeError(f"Unknown mesh_type={self.mesh_type}")
        
        mesh_width = mesh_nodes[1:] - mesh_nodes[:-1]

        assert len(mesh_nodes) == self.n_elem + 1
        assert len(mesh_width) == self.n_elem
        return mesh_nodes, mesh_width


    @staticmethod
    def generate_fe_nodes(boundaries_nodes: np.ndarray,
                          interp_nodes: np.ndarray) -> np.ndarray:
        """
        Map reference nodes from [-1,1] to each physical element [r_i, r_{i+1}],
        remove duplicate interface nodes, and return the global FE nodes.

        Parameters
        ----------
        boundaries_nodes : np.ndarray, shape (num_elements+1,)
            Coordinates of boundaries of finite elements (monotonically increasing or decreasing).
        interp_nodes : np.ndarray, shape (poly_order+1,)
            Interpolation nodes in [-1,1], e.g. Lobatto nodes.

        Returns
        -------
        global_nodes : np.ndarray, shape (num_elements * poly_order + 1,)
            Concatenated FE nodes over [0, domain_size] without duplicate interface nodes.
        """
        assert isinstance(boundaries_nodes, np.ndarray)
        assert isinstance(interp_nodes, np.ndarray)

        # reverse the nodes if the nodes are decreasing
        if np.all(np.diff(boundaries_nodes) < 0.0):
            boundaries_nodes = boundaries_nodes[::-1]
        if np.all(np.diff(interp_nodes) < 0.0):
            interp_nodes = interp_nodes[::-1]

        # check if the nodes are monotonically increasing
        assert np.all(np.diff(boundaries_nodes) > 0.0)
        assert np.all(np.diff(interp_nodes) > 0.0)


        # check if interp_nodes includes -1 and 1 at both boundaries
        assert np.isclose(interp_nodes[0], -1.0)
        assert np.isclose(interp_nodes[-1], 1.0)

        num_elements = boundaries_nodes.size - 1
        poly_order = interp_nodes.size - 1

        left = boundaries_nodes[:-1]
        right = boundaries_nodes[1:]

        # affine map
        mapped = 0.5*(right - left)[:, None] * interp_nodes[None, :] \
                 + 0.5*(right + left)[:, None]

        global_nodes = mapped.ravel()

        # remove duplicates: last node of each element except the last
        remove_idx = np.arange(poly_order,
                               (poly_order+1)*(num_elements-1)+poly_order,
                               poly_order+1)
        global_nodes = np.delete(global_nodes, remove_idx)

        # sanity check
        expected_len = num_elements * poly_order + 1
        assert global_nodes.size == expected_len

        return global_nodes


    @staticmethod
    def refine_interpolation_nodes(base_nodes: np.ndarray) -> np.ndarray:
        """
        Refine interpolation nodes by inserting midpoints between existing nodes
        and adding extra midpoints near the two boundaries.

        Parameters
        ----------
        base_nodes : np.ndarray
            Original interpolation nodes (e.g., Lobatto nodes on [-1,1]).

        Returns
        -------
        refined_nodes : np.ndarray
            Refined nodes, containing original nodes, all midpoints,
            and extra midpoints near the left and right boundaries.
            The output is returned in descending order (consistent with some FE codes).
        """
        assert isinstance(base_nodes, np.ndarray)
        assert base_nodes.ndim == 1
        
        # reverse the nodes if the nodes are decreasing
        monotonically_increasing_flag = True
        if np.all(np.diff(base_nodes) < 0.0):
            base_nodes = base_nodes[::-1]
            monotonically_increasing_flag = False
        assert np.all(np.diff(base_nodes) > 0.0)

        # insert midpoints between consecutive nodes
        midpoints = base_nodes[:-1] + 0.5*np.diff(base_nodes)
        refined_nodes = np.sort(np.concatenate((base_nodes, midpoints)))

        # extra midpoint near left boundary
        refined_nodes = np.insert(refined_nodes, 1, 0.5*(refined_nodes[0] + refined_nodes[1]))

        # # extra midpoint near right boundary
        # In the latest implementation, we don't need to add extra midpoint near right boundary
        # refined_nodes = np.insert(refined_nodes, len(refined_nodes)-1, 0.5*(refined_nodes[-1] + refined_nodes[-2]))

        # flip to descending order (optional, to match legacy FE codes)
        if monotonically_increasing_flag:
            refined_nodes = np.flip(refined_nodes)

        return refined_nodes


    @staticmethod
    def map_quadrature_to_physical_elements(
        boundaries_nodes : np.ndarray,
        interp_nodes     : np.ndarray,
        interp_weights   : np.ndarray,
        flatten          : bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Affine-map Gauss-Legendre (or any reference) quadrature nodes/weights from [-1,1]
        to each physical element [x_i, x_{i+1}] and return global points/weights.

        Parameters
        ----------
        boundaries_nodes : (N_e+1,) array_like
            Monotonically nondecreasing element boundaries [x_0, ..., x_{N_e}].
        interp_nodes : (q,) array_like
            Quadrature nodes in [-1,1] on the reference interval.
        interp_weights : (q,) array_like
            Corresponding quadrature weights on [-1,1].
        flatten : bool, default True
            If True, return 1D arrays of length N_e*q (global concatenation).
            If False, return 2D arrays of shape (N_e, q) per element.

        Returns
        -------
        quad_x : ndarray
            Mapped quadrature nodes. Shape (N_e*q,) if flatten else (N_e, q).
        quad_w : ndarray
            Mapped quadrature weights. Shape (N_e*q,) if flatten else (N_e, q).

        Notes
        -----
        Mapping formula for each element [xL, xR]:
            x = 0.5*(xR - xL)*xi + 0.5*(xR + xL)
            w = 0.5*(xR - xL)*w_ref
        """
        assert isinstance(boundaries_nodes, np.ndarray)
        assert isinstance(interp_nodes, np.ndarray)
        assert isinstance(interp_weights, np.ndarray)

        xb = np.asarray(boundaries_nodes, dtype=float).reshape(-1)
        xi = np.asarray(interp_nodes, dtype=float).reshape(-1)
        wi = np.asarray(interp_weights, dtype=float).reshape(-1)

        # basic checks
        assert xb.size >= 2
        assert xi.size == wi.size
        assert xi.size >= 1
        assert np.all(np.diff(xb) >= 0.0)


        Ne = xb.size - 1
        q = xi.size

        xL = xb[:-1]                    # (Ne,)
        xR = xb[1:]                     # (Ne,)
        h  = (xR - xL)[:, None]         # (Ne,1)

        # map nodes & weights: (Ne,q)
        quad_x = 0.5 * h * xi[None, :] + 0.5 * (xR + xL)[:, None]
        quad_w = 0.5 * h * wi[None, :]

        if flatten:
            return quad_x.reshape(-1), quad_w.reshape(-1)
        return quad_x, quad_w


    @staticmethod
    def fe_flat_to_block2d(
        flat   : np.ndarray,
        n_elem : int,
        endpoints_shared: bool,
    ) -> np.ndarray:

        r"""
        Convert a 1D concatenated FE grid into a 2D array of shape (n_elem, m),
        where m is the number of points per element.

            Two supported layouts for 'flat':
            1) endpoints_shared=True:
                Global array stores unique boundary nodes, so consecutive elements share 1 node.
                Expected length: len(flat) = n_elem * (m - 1) + m.
            2) endpoints_shared=False:
                Elements are simply stacked; no overlapping boundary nodes.
                Expected length: len(flat) = n_elem * m.

            Parameters
            ----------
            flat : (L,) 1D ndarray
                Concatenated points on the reference line.
            n_elem : int
                Number of finite elements.
            points_per_elem : Optional[int]
                Points per element (m). If None, it will be inferred from len(flat), n_elem,
                and endpoints_shared (if provided) or by trying both layouts.
            endpoints_shared : Optional[bool]
                If True/False, use the specified layout. If None, auto-detect by matching lengths.

            Returns
            -------
            blocks : (n_elem, m) ndarray
                The 2D array where each row corresponds to one element's local grid.

            Examples
            --------
            # Example 1: endpoints shared (unique global nodes)
            
                Suppose each element has m=4 nodes and n_elem=3:
                element 0: [0, 1, 2, 3]
                element 1: [3, 4, 5, 6]
                element 2: [6, 7, 8, 9]
                flat1: 
                    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]  (length = 3*4 - 2 = 10)
                blocks1:
                    [[0, 1, 2, 3],
                     [3, 4, 5, 6],
                     [6, 7, 8, 9]]

            # Example 2: endpoints NOT shared (stacked)
                Suppose each element has m=4 nodes and n_elem=3:
                element 0: [0, 1, 2, 3]
                element 1: [4, 5, 6, 7]
                element 2: [8, 9,10,11]
                flat2:
                    [0,1,2,3,  4,5,6,7,  8,9,10,11] (length = 3*4 = 12)
                blocks2:
                    [[0, 1, 2, 3],
                     [4, 5, 6, 7],
                     [8, 9,10,11]]
        """

        assert isinstance(flat, np.ndarray)
        assert flat.ndim == 1
        assert isinstance(n_elem, int)
        assert isinstance(endpoints_shared, bool)

        # ---- Construct the (n_elem, m) view/array ----
        if endpoints_shared:
            assert (flat.size - 1) % n_elem == 0
            # Overlapping windows with stride (m-1)
            m = (flat.size - 1) // n_elem + 1
            out = np.empty((n_elem, m), dtype=flat.dtype)
            for i_elem in range(n_elem):
                start = i_elem * (m - 1)
                out[i_elem, :] = flat[start : start + m]
            return out

        else:
            assert flat.size % n_elem == 0
            # Simple reshape when elements are stacked without overlap
            return flat.reshape(n_elem, -1)



class LagrangeShapeFunctions:
    """
    Lagrange shape functions (and their derivatives) built from a given set of nodes.
    This class replaces the ad-hoc broadcasting logic in
    - lagrange_polynomial(...)
    - lagrange_polynomial_t_obtain_quantities_uniform_grid(...)

    It supports:
      1) Reference-element evaluation: given reference nodes xj in [-1,1], evaluate Phi,dPhi at xi.
      2) Physical-element evaluation: given element-wise physical nodes rj (Ne,p+1) and eval points re (Ne,n_eval),
         evaluate Phi(re), dPhi(re) w.r.t. the physical coordinate r (chain rule handled).
      3) Interpolation to a global uniform grid with per-element selection (returns lag_poly and raw_indexes).
    """



    # @classmethod
    # def lagrange_basis_and_derivatives(cls, x_node: np.ndarray, x_eval: np.ndarray, atol: float = 1e-12) -> Tuple[np.ndarray, np.ndarray]:
    #     return cls.lagrange_basis_and_derivatives_old(x_node, x_eval, atol)


    @staticmethod
    def lagrange_basis_and_derivatives(x_node: np.ndarray, x_eval: np.ndarray, atol: float = 1e-20) -> Tuple[np.ndarray, np.ndarray]:
        r"""
        Compute Lagrange basis values and their first derivatives on each FE element,
        fully vectorized over elements, evaluation points, and basis indices.

        Parameters
        ----------
        x_node : ndarray, shape (n_elem, n_node)
            Element-local nodal coordinates x_j for each element (one row per element).
            For an element with polynomial degree p, m = p+1.
        x_eval : ndarray, shape (n_elem, n_eval)
            Evaluation points on each element (one row per element). 
            These can be quadrature points or the nodes themselves.
        atol : float
            Absolute tolerance to decide if an evaluation point exactly coincides
            with a node (used to avoid 0/0 and to apply the exact nodal limits).


        Returns
        -------
        L : ndarray, shape (n_elem, n_eval, n_node)
            Lagrange basis values. L[i, j, k] = L_k(x_eval_ij) on element i.
        dLdx : ndarray, shape (n_elem, n_eval, n_node)
            First derivatives of the Lagrange basis. dLdx[i, j, k] = d/dx L_k at x_eval_ij.


        Mathematics (per element)
        -------------------------
        Let nodes be {x_k}_{k=0}^{m-1}, m = n_node.
        Barycentric weights:
            ω_k = 1 / prod_{t≠k} (x_k - x_t)

        For a non-nodal evaluation point x:
            L_k(x) = (ω_k / (x - x_k)) / ( sum_t(ω_t / (x - x_t)) )

        Derivative at a non-nodal x:
            L'_k(x) = L_k(x) * ( sum_t(1/(x - x_t)) - 1/(x - x_k) )

        For nodal x = x_a, the derivative row equals the a-th row of the
        nodal differentiation matrix D, with off-diagonal entries
            D[a,b] = (c_a / c_b) / (x_a - x_b),  a ≠ b,
        where c_k = prod_{t≠k} (x_k - x_t), and diagonal fixed by
            D[a,a] = - sum_{b≠a}( D[a,b] ).


        Notes
        -----
        * Uses absolute matching tolerance `atol` with rtol=0 to detect x exactly at nodes.
        * Avoids forming large/small products at evaluation points (stable vs. "full-product then divide").
        """

        # basic checks
        assert isinstance(x_node, np.ndarray)
        assert isinstance(x_eval, np.ndarray)
        if x_node.ndim == 1:
            x_node = x_node[None, :]
        if x_eval.ndim == 1:
            x_eval = x_eval[None, :]
        assert x_node.ndim == 2
        assert x_eval.ndim == 2
        assert x_node.shape[0] == x_eval.shape[0]
        assert x_node.shape[1] >= 2
        assert x_eval.shape[1] >= 1

        # shapes
        n_elem, n_node = x_node.shape
        n_elem, n_eval = x_eval.shape

        # 1) Precompute barycentric weights: omega = 1 / c, where
        #    c_k = prod_{t≠k} (x_k - x_t), done per element.

        # diffs_nodes[i, k, t] = x_k - x_t (on element i)
        diffs_nodes = x_node[:, :, None] - x_node[:, None, :]                 # (n_elem, n_node, n_node)
        mask_offdiag = ~np.eye(n_node, dtype=bool)                            # (n_node, n_node) True for a≠b

        # c[i, k] = prod_{t≠k} (x_k - x_t): set diagonal to 1 so the product ignores it
        c = np.where(mask_offdiag[None, :, :], diffs_nodes, 1.0).prod(axis=2) # (n_elem, n_node)
        omega = 1.0 / c                                                       # (n_elem, n_node)

        # 2) Distances from evaluation points to nodes: dx = x - x_k
        dx = x_eval[:, :, None] - x_node[:, None, :]                          # (n_elem, n_eval, n_node)

        # 3) Identify nodal evaluations (x == x_k) using absolute tolerance
        is_nodal = np.isclose(x_eval[:, :, None], x_node[:, None, :], atol=atol, rtol=0.0)  # (n_elem, n_eval, n_node)
        has_nodal_row = is_nodal.any(axis=2)                                      # (n_elem, n_eval)

        # 4) Basis values at non-nodal points via barycentric formula
        #    L_k(x) = (omega_k /(x - x_k)) / sum_t (omega_t /(x - x_t))
        with np.errstate(divide='ignore', invalid='ignore'):
            inv_dx = 1.0 / dx                                               # (n_elem, n_point, n_node)
            numer = omega[:, None, :] * inv_dx                              # (n_elem, n_point, n_node)
            denom = np.sum(numer, axis=2, keepdims=True)                    # (n_elem, n_point, 1)
            L = numer / denom                                               # (n_elem, n_point, n_node)

        # Enforce nodal rows to be exact one-hot (interpolation property)
        # First zero out rows that contain nodal matches; then set matching entries to 1.
        L = np.where(np.isfinite(L), L, 0.0)
        L = np.where(has_nodal_row[:, :, None], 0.0, L)
        L = np.where(is_nodal, 1.0, L)

        # 5) Derivatives at non-nodal points:
        #    L'_k(x) = L_k(x) * ( sum_t 1/(x - x_t) - 1/(x - x_k) )
        with np.errstate(divide='ignore', invalid='ignore'):
            harmonic_sum = np.sum(inv_dx, axis=2, keepdims=True)            # (n_elem, n_point, 1)
            dLdx = L * (harmonic_sum - inv_dx)                              # (n_elem, n_point, n_node)

        # 6) Nodal differentiation matrix D on the nodes {x_k} (per element)
        #    Off-diagonal: D[a,b] = (c_a / c_b) / (x_a - x_b)
        #    Diagonal:     D[a,a] = -sum_{b≠a} D[a,b]
        with np.errstate(divide='ignore', invalid='ignore'):
            D = (c[:, :, None] / c[:, None, :]) / diffs_nodes               # (n_elem, n_node, n_node)
        D = np.where(mask_offdiag[None, :, :], D, 0.0)
        D[..., np.arange(n_node), np.arange(n_node)] = -np.sum(D, axis=2)

        # Overwrite derivative rows where x coincides with a node:
        # For each (i, j) with a nodal match, find which node index k* and set dLdx[i, j, :] = D[i, k*, :]
        if has_nodal_row.any():
            nodal_index = np.argmax(is_nodal, axis=2)                       # (n_elem, n_point), argmax is safe since rows with no nodal are masked by has_nodal_row
            i_idx, j_idx = np.nonzero(has_nodal_row)                        # indices of rows to overwrite
            k_idx = nodal_index[i_idx, j_idx]
            dLdx[i_idx, j_idx, :] = D[i_idx, k_idx, :]

        # 7) Cleanup any residual non-finite entries (should be rare)
        L = np.where(np.isfinite(L), L, 0.0)
        dLdx = np.where(np.isfinite(dLdx), dLdx, 0.0)

        return L, dLdx


    @staticmethod
    def _barycentric_weights(xj: np.ndarray) -> np.ndarray:
        """Compute barycentric weights w_j = 1 / ∏_{m≠j} (x_j - x_m)."""
        xj = np.asarray(xj, dtype=float)
        n = xj.size
        w = np.ones(n, dtype=float)
        for j in range(n):
            diff = xj[j] - xj
            diff[j] = 1.0  # avoid zero
            w[j] = 1.0 / np.prod(diff)
        return w




    @staticmethod
    def _derivative_matrix_at_nodes(xj: np.ndarray, w: np.ndarray) -> np.ndarray:
        """
        Derivatives of Lagrange basis evaluated at the nodes:
            D[k,j] = φ_j'(x_k)
        with φ_j the cardinal basis through xj. Closed form:
          for j≠k: D[k,j] = w[j] / ( w[k] * (x_k - x_j) )
          for j=k: D[k,k] = -∑_{m≠k} D[k,m]
        """
        xj = np.asarray(xj, dtype=float)
        n = xj.size
        D = np.empty((n, n), dtype=float)
        for k in range(n):
            for j in range(n):
                if j == k:
                    continue
                D[k, j] = w[j] / (w[k] * (xj[k] - xj[j]))
            D[k, k] = -np.sum(D[k, :]) + D[k, k] if np.isfinite(D[k, k]) else -np.sum(D[k, :])
        return D





@dataclass(frozen=True)
class GridData:
    """
    Immutable grid information needed for XC calculations.
    
    Parameters
    ----------
    finite_element_number : int
        Number of finite elements
    physical_nodes : np.ndarray
        Physical FE nodes
    quadrature_nodes : np.ndarray
        Quadrature points
    quadrature_weights : np.ndarray
        Quadrature weights
    """
    finite_element_number : int
    physical_nodes        : np.ndarray
    quadrature_nodes      : np.ndarray
    quadrature_weights    : np.ndarray
    
    # Deprecated: for backward compatibility
    @property
    def number_of_finite_elements(self) -> int:
        """Deprecated: use finite_element_number instead."""
        return self.finite_element_number


    def __post_init__(self):
        # check if the input parameters are valid
        assert self.finite_element_number > 0
        assert self.physical_nodes.ndim == 1
        assert self.quadrature_nodes.ndim == 1
        assert self.quadrature_weights.ndim == 1
        assert self.quadrature_nodes.shape[0] == self.quadrature_weights.shape[0]

    @classmethod
    def from_basic(
        cls, 
        domain_size             : float,
        finite_element_number   : int,
        polynomial_order        : int,
        quadrature_point_number : int,
        mesh_type               : str,
        mesh_concentration      : float,
    ) -> 'GridData':

        # Generate Lobatto interpolation nodes on reference interval [-1, 1]
        interp_nodes_ref, _ = Quadrature1D.lobatto(polynomial_order)
        
        # Generate mesh boundaries
        mesh1d = Mesh1D(
            domain_size            = domain_size,
            finite_elements_number = finite_element_number,
            mesh_type              = mesh_type,
            clustering_param       = mesh_concentration,
        )
        boundaries_nodes, _ = mesh1d.generate_mesh_nodes_and_width()

        # Generate standard FE nodes
        global_nodes = Mesh1D.generate_fe_nodes(
            boundaries_nodes = boundaries_nodes,
            interp_nodes     = interp_nodes_ref
        )

        # Generate Gauss-Legendre quadrature nodes and weights
        quadrature_nodes_ref, quadrature_weights_ref = Quadrature1D.gauss_legendre(quadrature_point_number)

        # Map quadrature to physical elements
        quadrature_nodes, quadrature_weights = Mesh1D.map_quadrature_to_physical_elements(
            boundaries_nodes = boundaries_nodes,
            interp_nodes     = quadrature_nodes_ref,
            interp_weights   = quadrature_weights_ref,
            flatten          = True
        )

        return cls(
            finite_element_number = finite_element_number,
            physical_nodes        = global_nodes,
            quadrature_nodes      = quadrature_nodes,
            quadrature_weights    = quadrature_weights,
        )
        



class RadialOperatorsBuilder:
    """
    Assemble radial FE operators and interpolation matrices.

    Initialize from individual parameters, or use from_grid_data() for GridData.
    """

    def __init__(self, 
        finite_element_number     : Optional[int]        = None,
        physical_nodes            : Optional[np.ndarray] = None, 
        quadrature_nodes          : Optional[np.ndarray] = None,
        quadrature_weights        : Optional[np.ndarray] = None,
        verbose                   : bool = False,
        builder_label             : str = "",  # label for the builder, will be printed in the summary

        # Deprecated parameters
        number_of_finite_elements : Optional[int] = None,  # Deprecated: use finite_element_number instead
    ):
        """
        Initialize radial operators builder from individual parameters.
        
        For initialization from GridData, use RadialOperatorsBuilder.from_grid_data().
        
        Parameters
        ----------
        finite_element_number : int
            Number of finite elements
        physical_nodes : np.ndarray
            Physical FE nodes
        quadrature_nodes : np.ndarray
            Quadrature points
        quadrature_weights : np.ndarray
            Quadrature weights
        verbose : bool, optional
            If True, print initialization summary. Default: False
        builder_label : str, optional
            Label for the builder. Default: ""
        number_of_finite_elements : int, optional, deprecated
            Deprecated parameter. Use 'finite_element_number' instead.
        """

        # Handle deprecated parameter
        if number_of_finite_elements is not None:
            if finite_element_number is not None:
                raise RuntimeError("Cannot specify both 'finite_element_number' and deprecated 'number_of_finite_elements'. Use 'finite_element_number' only.")
            finite_element_number = number_of_finite_elements
            if verbose:
                print("NUMBER_OF_FINITE_ELEMENTS_DEPRECATED_WARNING")

        # Required parameters (must be at the beginning)
        assert physical_nodes is not None
        assert quadrature_nodes is not None
        assert quadrature_weights is not None
        assert finite_element_number is not None
        
        # dimension and shape checks
        assert physical_nodes.ndim == 1
        assert quadrature_nodes.ndim == 1
        assert quadrature_weights.ndim == 1
        assert quadrature_nodes.shape[0] == quadrature_weights.shape[0]

        self.assembly_dtype        = MATRIX_ASSEMBLY_DTYPE
        self.finite_element_number = finite_element_number
        self.physical_nodes        = np.asarray(physical_nodes)
        self.quadrature_nodes      = np.asarray(quadrature_nodes)
        self.quadrature_weights    = np.asarray(quadrature_weights)
        self.verbose               = verbose
        
        # Reshape to element-wise structure
        self._reshape_grid_data()
        
        # Compute Lagrange basis functions
        self._compute_basis_functions()
        
        # Print summary if verbose
        if self.verbose:
            self._print_initialization_summary(builder_label)
    
    
    @classmethod
    def from_grid_data(cls, grid_data: GridData, verbose: bool = False, builder_label: str = "") -> 'RadialOperatorsBuilder':
        """
        Create RadialOperatorsBuilder from GridData.
        
        Convenience factory method for cleaner code.
        
        Parameters
        ----------
        grid_data : GridData
            Grid data object
        verbose : bool, optional
            If True, print initialization summary. Default: False
        builder_label : str, optional
            Label for the builder. Default: ""
        
        Returns
        -------
        RadialOperatorsBuilder
            Initialized operators builder
        
        Example
        -------
        >>> grid_data = GridData(...)
        >>> ops = RadialOperatorsBuilder.from_grid_data(grid_data, verbose=True)
        """
        return cls(
            finite_element_number = grid_data.finite_element_number,
            physical_nodes        = grid_data.physical_nodes,
            quadrature_nodes      = grid_data.quadrature_nodes,
            quadrature_weights    = grid_data.quadrature_weights,
            verbose               = verbose,
            builder_label         = builder_label,
        )


    def _reshape_grid_data(self):
        """Reshape 1D arrays to (n_elem, n_points) structure."""
        self.physical_nodes_reshaped = Mesh1D.fe_flat_to_block2d(
            self.physical_nodes, 
            self.finite_element_number, 
            endpoints_shared=True
        )
        self.quadrature_nodes_reshaped = Mesh1D.fe_flat_to_block2d(
            self.quadrature_nodes, 
            self.finite_element_number, 
            endpoints_shared=False
        )
        self.quadrature_weights_reshaped = Mesh1D.fe_flat_to_block2d(
            self.quadrature_weights, 
            self.finite_element_number, 
            endpoints_shared=False
        )
        
        self.physical_node_number = self.physical_nodes_reshaped.shape[1]
        self.quadrature_node_number = self.quadrature_nodes_reshaped.shape[1]
    
    
    def _compute_basis_functions(self):
        """Compute Lagrange basis functions and derivatives."""
        self.lagrange_basis, self.lagrange_basis_derivatives = \
            LagrangeShapeFunctions.lagrange_basis_and_derivatives(
                x_node=self.physical_nodes_reshaped,
                x_eval=self.quadrature_nodes_reshaped
            )
        # Shape: (n_elem, n_quad, n_basis)
    
    
    def _print_initialization_summary(self, label: str = ""):
        """Print initialization summary."""
        print("=" * 75)
        print("RadialOperatorsBuilder ({})".format(label).center(75))
        print("=" * 75)
        print(f"\t Number of elements            : {self.finite_element_number}")
        print(f"\t Physical nodes per element    : {self.physical_nodes_reshaped.shape[1]}")
        print(f"\t Quadrature points per element : {self.quadrature_node_number}")
        print(f"\t Total DOFs                    : {self.finite_element_number * self.physical_node_number + 1}")
        print(f"\t Total quadrature points       : {len(self.quadrature_nodes)}")
        print(f"\t Lagrange basis shape          : {self.lagrange_basis.shape}")
        print(f"\t Lagrange derivatives shape    : {self.lagrange_basis_derivatives.shape}")
        print()
    
    # Deprecated: for backward compatibility
    @property
    def number_of_finite_elements(self) -> int:
        """Deprecated: use finite_element_number instead."""
        return self.finite_element_number
    
    
    def get_H_kinetic(self) -> np.ndarray:
        """
        Kinetic energy matrix
        """
        # check if the kinetic energy matrix has been computed
        if hasattr(self, "_H_kinetic"):
            return self._H_kinetic
        
        # compute the kinetic energy matrix
        H_kinetic = - 0.5 * self.get_laplacian()
        
        # store the kinetic energy matrix
        self._H_kinetic = H_kinetic

        return H_kinetic


    def get_laplacian(self) -> np.ndarray:
        """
        Laplacian matrix
        """
        # check if the Laplacian matrix has been computed
        if hasattr(self, "_laplacian"):
            return self._laplacian
        
        # compute the Laplacian matrix
        laplacian_local = - np.einsum("emi,emk,em->eik", 
            self.lagrange_basis_derivatives, 
            self.lagrange_basis_derivatives, 
            self.quadrature_weights_reshaped, 
            optimize=True)

        # Assemble local matrices into global matrix
        laplacian = self._assemble_local_to_global_matrix(laplacian_local)
        
        # store the Laplacian matrix
        self._laplacian = laplacian
        
        return laplacian


    def build_potential_matrix(self, potential_values: np.ndarray) -> np.ndarray:
        """
        Construct Hamiltonian matrix from potential energy values at quadrature points.
        
        This is a general-purpose matrix builder that can be used for any potential:
        - Nuclear Coulomb potential: V(r) = -Z/r
        - Local pseudopotential: V_loc(r)
        - Hartree potential: V_H(r)
        - Exchange-correlation potential: V_xc(r)
        
        Theory
        ------
        Given potential V(r) sampled at quadrature points, compute:
            H_V[i,j] = ∫ φ_i(r) V(r) φ_j(r) dr
        
        Using finite element basis φ and Gauss quadrature.
        
        Parameters
        ----------
        potential_values : np.ndarray
            Potential energy values at quadrature points.
            Can be 1D (flat, length = n_elem * n_quad) or 
            2D (reshaped, shape = (n_elem, n_quad))
        
        Returns
        -------
        H_potential : np.ndarray
            Assembled global Hamiltonian matrix contribution from this potential.
            Shape: (n_global_dofs, n_global_dofs)
        
        Examples
        --------
        >>> # Nuclear Coulomb potential
        >>> V_nuc = -Z / ops.quadrature_nodes
        >>> H_nuc = ops.build_potential_matrix(V_nuc)
        
        >>> # XC potential
        >>> V_xc = compute_xc_potential(rho)
        >>> H_xc = ops.build_potential_matrix(V_xc)
        """
        # Reshape potential if needed
        if potential_values.ndim == 1:
            assert potential_values.size == self.quadrature_nodes_reshaped.size
            potential_reshaped = potential_values.reshape(self.quadrature_nodes_reshaped.shape)
        elif potential_values.ndim == 2:
            assert potential_values.shape[0] == self.finite_element_number
            assert potential_values.shape[1] == self.quadrature_node_number
            potential_reshaped = potential_values
        else:
            raise RuntimeError("POTENTIAL_VALUES_NDIM_ERROR")
        
        # Compute local element matrices: H^e[i,k] = ∫ φ_i V φ_k dr
        H_local = np.einsum("emi,emk,em,em->eik", 
            self.lagrange_basis,       # φ_i at quadrature points
            self.lagrange_basis,       # φ_k at quadrature points
            self.quadrature_weights_reshaped,  # quadrature weights
            potential_reshaped,        # V(r) at quadrature points
            optimize=True)
        
        # Assemble local matrices into global matrix
        H_potential = self._assemble_local_to_global_matrix(H_local)
        
        return H_potential


    def build_metagga_kinetic_density_matrix(self, de_xc_dtau: np.ndarray) -> np.ndarray:
        """
        Construct Hamiltonian matrix from meta-GGA kinetic density term at quadrature points.
        
        This implements the meta-GGA contribution to the Hamiltonian:
            H_metagga[i,j] = ∫ φ'_i * (0.5 * de_xc_dtau) * φ'_j dr
                          + ∫ φ_i * (0.5 * w/r * d(de_xc_dtau)/dr) * φ_j dr
        
        Theory
        ------
        For meta-GGA functionals, the kinetic energy density τ contributes additional
        terms to the Hamiltonian. The implementation follows:
            V3 = de_xc_dtau
            V3grad = d(V3)/dr (computed via derivative matrix)
            
            Term 1: ∫ φ'_i * (0.5 * V3 * w) * φ'_j dr
            Term 2: ∫ φ_i * (0.5 * w/r * V3grad) * φ_j dr
        
        Parameters
        ----------
        de_xc_dtau : np.ndarray
            Derivative of XC energy density w.r.t. kinetic energy density τ.
            Can be 1D (flat, length = n_elem * n_quad) or 
            2D (reshaped, shape = (n_elem, n_quad))
        
        Returns
        -------
        H_metagga : np.ndarray
            Assembled global Hamiltonian matrix contribution from meta-GGA kinetic density term.
            Shape: (n_global_dofs, n_global_dofs)
        """
        # Reshape de_xc_dtau if needed
        if de_xc_dtau.ndim == 1:
            assert de_xc_dtau.size == self.quadrature_nodes_reshaped.size
            de_xc_dtau_reshaped = de_xc_dtau.reshape(self.quadrature_nodes_reshaped.shape)
        elif de_xc_dtau.ndim == 2:
            assert de_xc_dtau.shape[0] == self.finite_element_number
            assert de_xc_dtau.shape[1] == self.quadrature_node_number
            de_xc_dtau_reshaped = de_xc_dtau
        else:
            raise RuntimeError("DE_XC_DTAU_NDIM_ERROR")
        
        # Get derivative matrix for computing d(de_xc_dtau)/dr
        D = self.get_derivative_matrix_with_quadrature_basis()  # Shape: (n_elem, n_quad, n_quad)
        
        # Compute gradient: v3grad = D @ de_xc_dtau (element-wise)
        # For each element: v3grad[e, i] = sum_k D[e, i, k] * de_xc_dtau[e, k]
        de_xc_dtau_grad = np.einsum("eik,ek->ei", D, de_xc_dtau_reshaped)
        
        # Term 1: ∫ φ'_i * (0.5 * de_xc_dtau * w) * φ'_j dr
        # Using derivative basis functions
        H_term1_local = np.einsum("emi,emk,em,em->eik",
            self.lagrange_basis_derivatives,  # φ'_i
            self.lagrange_basis_derivatives,  # φ'_j
            self.quadrature_weights_reshaped,  # w
            0.5 * de_xc_dtau_reshaped,  # 0.5 * de_xc_dtau
            optimize=True)
        
        # Term 2: ∫ φ_i * (0.5 * w/r * d(de_xc_dtau)/dr) * φ_j dr
        # Using regular basis functions
        # Note: 0.5 * w / r * v3grad
        weight_over_r = 0.5 * self.quadrature_weights_reshaped / self.quadrature_nodes_reshaped
        H_term2_local = np.einsum("emi,emk,em,em->eik",
            self.lagrange_basis,  # φ_i
            self.lagrange_basis,  # φ_j
            weight_over_r,  # 0.5 * w / r
            de_xc_dtau_grad,  # d(de_xc_dtau)/dr
            optimize=True)
        
        # Combine both terms
        H_metagga_local = H_term1_local + H_term2_local
        
        # Assemble local matrices into global matrix
        H_metagga = self._assemble_local_to_global_matrix(H_metagga_local)
        
        return H_metagga




    def get_nuclear_coulomb_potential(self, z_nuclear: float) -> np.ndarray:
        """
        Compute nuclear Coulomb potential: V_nuc(r) = -Z/r
        
        Parameters
        ----------
        z_nuclear : float
            Nuclear charge (atomic number for all-electron, or effective charge)
        
        Returns
        -------
        V_nuclear : np.ndarray
            Nuclear Coulomb potential at quadrature points.
            Shape: (n_elem * n_quad,) flat array
        
        Notes
        -----
        For all-electron calculations, z_nuclear is the atomic number Z.
        The negative sign accounts for attractive electron-nucleus interaction.
        """
        try:
            z_nuclear = float(z_nuclear)
        except:
            raise RuntimeError("Z_NUCLEAR_NOT_FLOAT_ERROR")
        
        # V_nuc(r) = -Z/r at all quadrature points
        V_nuclear = -z_nuclear / self.quadrature_nodes
        
        return V_nuclear


    def get_H_r_inv_sq(self) -> np.ndarray:
        """
        Inverse square of the radial distance matrix
        """
        # check if the inverse square of the radial distance matrix has been computed
        if hasattr(self, "_H_r_inv_sq"):
            return self._H_r_inv_sq
    
        # compute the inverse square of the radial distance matrix
        H_r_inv_sq_local = np.einsum("emi,emk,em,em->eik", 
            self.lagrange_basis, 
            self.lagrange_basis, 
            self.quadrature_weights_reshaped, 
            1.0 / self.quadrature_nodes_reshaped**2, 
            optimize=True)
        
        # Assemble local matrices into global matrix
        H_r_inv_sq = self._assemble_local_to_global_matrix(H_r_inv_sq_local)
        
        # store the inverse square of the radial distance matrix
        self._H_r_inv_sq = H_r_inv_sq
        
        return H_r_inv_sq


    def get_S(self, exclude_boundary: bool = False) -> np.ndarray:
        """
        Overlap matrix S
        """
        # check if the overlap matrix has been computed
        if hasattr(self, "_S"):
            if exclude_boundary:
                return self._S[1:-1,1:-1]
            return self._S
        
        # compute the overlap matrix
        S_local = np.einsum("emi,emk,em->eik", 
            self.lagrange_basis, 
            self.lagrange_basis, 
            self.quadrature_weights_reshaped, 
            optimize=True)

        # Assemble local matrices into global matrix
        S = self._assemble_local_to_global_matrix(S_local)
        
        # store the overlap matrix
        self._S = S

        if exclude_boundary:
            return S[1:-1,1:-1]

        return S


    def get_S_inv_sqrt(self) -> np.ndarray:
        """
        Inverse square root of the overlap matrix: S^(-1/2)
        Computed via eigendecomposition: S^(-1/2) = V * Λ^(-1/2) * V^T
        """
        if hasattr(self, "_S_inv_sqrt"):
            return self._S_inv_sqrt
        
        S = self.get_S()
        try:
            eigvals, eigvecs = np.linalg.eigh(S, UPLO='L')
        except TypeError:
            S64 = np.asarray(S, dtype=np.float64)
            eigvals, eigvecs = np.linalg.eigh(S64, UPLO='L')
        S_inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T
        S_inv_sqrt = 0.5 * (S_inv_sqrt + S_inv_sqrt.T)  # Symmetrize
        
        self._S_inv_sqrt = S_inv_sqrt
        return S_inv_sqrt


    def get_derivative_matrix(self) -> np.ndarray:
        """
        Differentiation matrix for computing derivatives at quadrature points.
        
        .. deprecated:: 
            This method is deprecated and will be removed in a future version.
            Use :meth:`get_derivative_matrix_with_quadrature_basis` instead.
        
        This matrix D enables direct computation of derivatives from function values:
            f'(x_quad) = D @ f(x_quad)
        
        Theory:
        -------
        For a function f(x) = sum_j c_j φ_j(x) in the FE space,
        we have f' = (dL/dx) @ L⁺ @ f, where:
        - L[i,j] = φ_j(x_quad_i) : basis functions at quadrature points
        - L⁺ : pseudoinverse of L (maps function values -> coefficients)
        - dL/dx[i,j] = φ'_j(x_quad_i) : basis derivatives at quadrature points
        
        Returns
        -------
        np.ndarray
            Shape: (n_elements, n_quad_points, n_quad_points)
            D[e, i, k] maps f(x_k) → f'(x_i) within element e
        """
        warnings.warn(
            "'get_derivative_matrix()' is deprecated and will be removed in a future version. "
            "Use 'get_derivative_matrix_with_quadrature_basis()' instead. "
            "The new method uses quadrature points as basis nodes, which is more efficient "
            "and matches the reference implementation.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # ============================================================================
        # OLD IMPLEMENTATION (preserved for reference, not executed)
        # ============================================================================
        # The original implementation using pseudoinverse of Lagrange basis.
        # This approach computes D = (dL/dx) @ L⁺ where L⁺ is the pseudoinverse.
        # Kept here for reference in case needed for future use cases.
        #
        # if hasattr(self, "_derivative_matrix"):
        #     return self._derivative_matrix
        # 
        # # Compute pseudoinverse of Lagrange basis: L⁺ maps values → coefficients
        # # Shape: (n_elem, n_basis, n_quad) where n_basis = physical_node_number
        # n_elem = self.number_of_finite_elements
        # n_basis = self.physical_node_number
        # n_quad = self.quadrature_node_number
        # 
        # basis_pseudoinverse = np.zeros((n_elem, n_basis, n_quad))
        # for elem_idx in range(n_elem):
        #     # self.lagrange_basis[elem_idx] has shape (n_quad, n_basis)
        #     basis_pseudoinverse[elem_idx] = np.linalg.pinv(self.lagrange_basis[elem_idx])
        # 
        # # Differentiation matrix: D = (dL/dx) @ L⁺
        # # Maps function values at quadrature points to derivative values
        # derivative_matrix = np.matmul(self.lagrange_basis_derivatives, basis_pseudoinverse)
        # 
        # self._derivative_matrix = derivative_matrix
        # return derivative_matrix
        
        # For backward compatibility, delegate to the new method
        # This ensures the function still works but uses the improved implementation
        return self.get_derivative_matrix_with_quadrature_basis()



    def get_derivative_matrix_with_quadrature_basis(self) -> np.ndarray:
        """
        Compute derivative matrix where basis function nodes = quadrature points.
        
        This method computes a special derivative matrix where the basis function nodes
        are set to be the same as quadrature points. This matches the reference code:
            D = lagrange_polynomial(r_quad1, q, r_quad1, q, S)[1]
        
        In this case:
        - n_quad == n_basis (quadrature points = basis nodes)
        - D[i, j] directly gives derivative of basis function j at quadrature point i
        - No need for inverse or pseudoinverse!
        
        The derivative matrix D satisfies:
            f'(x_quad) = D @ f(x_quad)
        
        where:
        - f(x_quad): function values at quadrature points, shape (n_elem, n_quad)
        - f'(x_quad): derivative values at quadrature points, shape (n_elem, n_quad)
        - D: derivative matrix, shape (n_elem, n_quad, n_quad)
        
        Returns
        -------
        np.ndarray
            Shape: (n_elements, n_quad_points, n_quad_points)
            D[e, i, k] maps f(x_k) → f'(x_i) within element e
        """
        if hasattr(self, "_derivative_matrix_quad_basis"):
            return self._derivative_matrix_quad_basis
        
        # Compute Lagrange basis and derivatives where nodes = quadrature points
        _, lagrange_basis_derivatives_quad = \
            LagrangeShapeFunctions.lagrange_basis_and_derivatives(
                x_node=self.quadrature_nodes_reshaped,  # Basis nodes = quadrature points
                x_eval=self.quadrature_nodes_reshaped,   # Evaluate at quadrature points
            )
        
        derivative_matrix = lagrange_basis_derivatives_quad  # (n_elem, n_quad, n_quad)
        
        self._derivative_matrix_with_quadrature_basis = derivative_matrix

        return derivative_matrix


    def get_derivative_matrix_with_physical_basis(self) -> np.ndarray:
        """
        Compute derivative matrix where basis function nodes = physical nodes.
        """
        if hasattr(self, "_derivative_matrix_with_physical_basis"):
            return self._derivative_matrix_with_physical_basis
        
        _, lagrange_basis_derivatives_physical = \
            LagrangeShapeFunctions.lagrange_basis_and_derivatives(
                x_node=self.physical_nodes_reshaped,
                x_eval=self.physical_nodes_reshaped,
            )
        
        # Store the derivative matrix with physical basis
        self._derivative_matrix_with_physical_basis = lagrange_basis_derivatives_physical
        return lagrange_basis_derivatives_physical



    def get_global_interpolation_matrix(self) -> np.ndarray:
        """
        Global interpolation matrix for evaluating functions at quadrature points.
        
        This matrix maps global nodal coefficients to function values at all quadrature points:
            f(x_quad) = interpolation_matrix @ f_nodal
        
        Theory
        ------
        For a function f(x) represented by nodal values on the global FE grid,
        this matrix evaluates f at all quadrature points across all elements.
        
        Unlike _assemble_local_to_global_matrix (which accumulates contributions),
        this constructs a rectangular interpolation matrix by:
        1. Placing local Lagrange basis matrices in block-diagonal form
        2. Merging columns corresponding to shared boundary nodes
        3. Removing duplicate columns
        
        Returns
        -------
        np.ndarray
            Shape: (n_elements * n_quad_points, n_global_nodes)
            where n_global_nodes = n_elements * polynomial_order + 1
            
            Example: n_elem=17, p=31, n_quad=95 → shape (1615, 528)
        
        Usage
        -----
        Essential for Hartree-Fock exchange (HF, PBE0) where orbital overlaps
        must be computed at quadrature points:
            psi_i(x_quad) = interp_matrix @ psi_i_nodal
        
        Not needed for LDA, GGA, or meta-GGA functionals.
        """
        if hasattr(self, "_global_interpolation_matrix"):
            return self._global_interpolation_matrix
        
        n_elem   = self.finite_element_number
        n_quad   = self.quadrature_node_number
        n_local  = self.physical_node_number
        n_global = n_elem * (self.physical_node_number - 1) + 1
        
        # Method: Build block-diagonal then remove duplicate columns
        # Alternative could use fancy indexing, but this is clearer
        
        # Reshape lagrange_basis from (n_elem, n_quad, n_local) to 2D blocks
        # Result shape: (n_elem * n_quad, n_elem * n_local)
        interp_matrix = np.zeros((n_elem * n_quad, n_elem * n_local))
        
        for elem_idx in range(n_elem):
            row_slice = slice(elem_idx * n_quad, (elem_idx + 1) * n_quad)
            col_slice = slice(elem_idx * n_local, (elem_idx + 1) * n_local)
            interp_matrix[row_slice, col_slice] = self.lagrange_basis[elem_idx]
        
        # Shared boundary nodes: right edge of each element (except last)
        # These correspond to duplicate columns that must be merged
        shared_cols = np.array([n_local * (i + 1) - 1 for i in range(n_elem - 1)])
        
        # Merge: add right neighbor column to left neighbor column
        interp_matrix[:, shared_cols] += interp_matrix[:, shared_cols + 1]
        
        # Remove duplicate columns
        interp_matrix = np.delete(interp_matrix, shared_cols + 1, axis=1)
        
        assert interp_matrix.shape == (n_elem * n_quad, n_global), \
            f"Shape mismatch: expected {(n_elem * n_quad, n_global)}, get {interp_matrix.shape}."
        
        self._global_interpolation_matrix = interp_matrix
        return interp_matrix


    def get_lagrange_basis_pseudoinverse(self) -> np.ndarray:
        """
        Compute the pseudoinverse of the Lagrange basis for each finite element.
        """
        if hasattr(self, "_lagrange_basis_pseudoinverse"):
            return self._lagrange_basis_pseudoinverse
        
        lagrange_basis_pseudoinverse = np.zeros(
            (self.finite_element_number, self.physical_node_number, self.quadrature_node_number)
        )
        for elem_idx in range(self.finite_element_number):
            lagrange_basis_pseudoinverse[elem_idx, :, :] = np.linalg.pinv(self.lagrange_basis[elem_idx])

        self._lagrange_basis_pseudoinverse = lagrange_basis_pseudoinverse

        return lagrange_basis_pseudoinverse



    def assemble_poisson_rhs_vector_no_bc(self, rho: np.ndarray) -> np.ndarray:
        """
        Assemble the right-hand side vector for the Poisson equation.
        
        Computes the RHS vector for solving d²u/dr² = -4πrρ(r) where u = rV.
        The weak form involves integrating -4πrρ(r) against test functions.
        
        Parameters
        ----------
        rho : np.ndarray
            Electron density at quadrature points
            Shape: (n_elements * n_quad_points,)
        
        Returns
        -------
        rhs_vector : np.ndarray
            Assembled global RHS vector
            Shape: (n_global_dofs,)
            Note: Caller must set boundary values:
                rhs_vector[0] = left_bc
                rhs_vector[-1] = right_bc
        
        Notes
        -----
        - The local RHS contribution for each element is:
                RHS_local[i] = ∫ φ_i(r) * (-4πrρ(r)) dr
            These are assembled into a global vector, with contributions from
            adjacent elements added at shared boundary nodes.
        """
        # Type and shape validation
        assert isinstance(rho, np.ndarray)
        assert rho.ndim == 1
        assert rho.shape[0] == self.finite_element_number * self.quadrature_node_number
        # Reshape density to element-wise structure
        rho_reshaped = rho.reshape(self.finite_element_number, self.quadrature_node_number)
        
        # Compute source term at quadrature points: -4πrρ
        r_rho = self.quadrature_nodes_reshaped * rho_reshaped
        source = - 4.0 * np.pi * r_rho

        
        # Compute local RHS: ∫ φ_i(r) * source(r) dr
        rhs_vector_local = np.einsum(
            "emi,em,em->ei",
            self.lagrange_basis,              # φ_i at quadrature points
            self.quadrature_weights_reshaped, # quadrature weights
            source,                           # -4πrρ at quadrature points
            optimize=True
        )  # Shape: (n_elem, n_physical_nodes)

        # Assemble local vectors into global vector
        # Shared boundary nodes will have contributions from adjacent elements added
        rhs_vector = self._assemble_local_to_global_vector(rhs_vector_local)

        return rhs_vector




    def assemble_laplacian_rhs_vector_with_gradient_and_no_bc(self, grad_rho: np.ndarray) -> np.ndarray:
        """
        Assemble the right-hand side vector for the Laplacian weak form.
        
        Weak form (integration by parts on [0, R]):
            ∫₀^R φ_i ℓ r² dr = -∫₀^R φ_i' r² ρ' dr + [φ_i(r) r² ρ'(r)]₀^R
        where ℓ = ∇²ρ (Laplacian of density), φ_i are test functions, ρ' = dρ/dr.
        
        This method implements only the first term on the RHS (the integral with
        negative sign): -∫₀^R φ_i' r² ρ' dr. The boundary term is handled elsewhere.
        
        Parameters
        ----------
        grad_rho : np.ndarray
            Radial derivative of density ρ' = dρ/dr at quadrature points.
        """
        # Type and shape validation
        assert isinstance(grad_rho, np.ndarray)
        assert grad_rho.ndim == 1
        assert grad_rho.shape[0] == self.finite_element_number * self.quadrature_node_number
        
        # Reshape gradient of density to element-wise structure
        grad_rho_reshaped = grad_rho.reshape(self.finite_element_number, self.quadrature_node_number)

        # Compute source term at quadrature points: - r² ρ'
        source = - self.quadrature_nodes_reshaped ** 2 * grad_rho_reshaped

        # Compute local RHS: ∫ φ_i(r) * source(r) dr
        rhs_vector_local = np.einsum(
            "emi,em,em->ei",
            self.lagrange_basis_derivatives,  # φ_i' at quadrature points
            self.quadrature_weights_reshaped, # quadrature weights
            source,                           # - r² ρ' at quadrature points
            optimize=True
        )  # Shape: (n_elem, n_physical_nodes)

        # Assemble local vectors into global vector
        # Shared boundary nodes will have contributions from adjacent elements added
        rhs_vector = self._assemble_local_to_global_vector(rhs_vector_local)

        return rhs_vector


    def assemble_laplacian_rhs_vector_with_rho_and_no_bc(self, rho: np.ndarray) -> np.ndarray:
        """
        Assemble the right-hand side vector for the Laplacian weaker form.

        Weaker form (integration by parts on [0, R]):
            ∫₀^R φ_i ℓ r² dr = ∫₀^R (r² φ_i')' ρ dr - [r² φ_i'(r) ρ(r)]₀^R + [r² φ_i(r) ρ'(r)]₀^R
        where ℓ = ∇²ρ (Laplacian of density), φ_i are test functions, ρ' = dρ/dr.

        Uses (r² φ_i')' = 2r φ_i' + r² φ_i'' at quadrature points.

        Parameters
        ----------
        rho : np.ndarray
            Electron density at quadrature points.
        grad_rho : np.ndarray
            Radial derivative of density ρ' = dρ/dr at quadrature points.
        """
        # check if the input parameters are valid
        assert isinstance(rho, np.ndarray)
        assert rho.ndim == 1
        assert rho.shape[0] == self.finite_element_number * self.quadrature_node_number

        # Reshape density to element-wise structure
        rho_reshaped = rho.reshape(self.finite_element_number, self.quadrature_node_number)

        # Compute local RHS: ∫ φ_i(r) * source(r) dr
        rhs_vector_local = np.einsum(
            "em, eml, el, eli -> ei",
            self.quadrature_weights_reshaped * rho_reshaped,  # (n_elements, n_quad)
            self.derivative_matrix,                           # (n_elements, n_quad, n_quad)
            self.quadrature_nodes_reshaped ** 2,              # (n_elements, n_quad)
            self.lagrange_basis_derivatives,                  # (n_elements, n_quad, n_physical_nodes)
            optimize=True
        )

        rhs_vector = self._assemble_local_to_global_vector(rhs_vector_local)

        return rhs_vector


    def _assemble_local_to_global_vector(self, local_vector: np.ndarray) -> np.ndarray:
        """
        Assemble local element vectors into a global vector.
        
        This method takes local element-wise vectors (shape: [n_elements, n_dofs_per_elem])
        and assembles them into a single global vector by adding overlapping contributions 
        at shared nodes.
        
        Parameters
        ----------
        local_vector : np.ndarray
            Local element vectors with shape (n_elements, n_local_dofs)
            where n_local_dofs = polynomial_order + 1 = physical_node_number
            
        Returns
        -------
        np.ndarray
            Assembled global vector with shape (n_global_dofs,)
            where n_global_dofs = n_elements * polynomial_order + 1
            
        Notes
        -----
        The assembly process accounts for shared endpoints between adjacent elements,
        accumulating contributions from all elements that share a degree of freedom.
        
        For example, with 3 elements and polynomial order 2:
            Element 0: [a0, a1, a2]  → global indices [0, 1, 2]
            Element 1: [b0, b1, b2]  → global indices [2, 3, 4]  (b0 adds to global[2])
            Element 2: [c0, c1, c2]  → global indices [4, 5, 6]  (c0 adds to global[4])
        
        Global vector: [a0, a1, a2+b0, b1, b2+c0, c1, c2]
        """
        # Initialize global vector
        global_size = self.finite_element_number * (self.physical_node_number - 1) + 1
        global_vector = np.zeros(global_size, dtype=self.assembly_dtype)
        
        # Get assembly indices (mapping from local to global DOFs)
        indices_global = self._build_assembly_indices_vector()
        
        # Assemble: add local contributions to global vector
        # np.add.at handles accumulation at repeated indices automatically
        np.add.at(global_vector, indices_global, np.asarray(local_vector, dtype=self.assembly_dtype).reshape(-1))
        
        return global_vector


    def _assemble_local_to_global_matrix(self, local_matrix: np.ndarray) -> np.ndarray:
        """
        Assemble local element matrices into a global matrix.
        
        This method takes local element-wise matrices (shape: [n_elements, n_dofs_per_elem, n_dofs_per_elem])
        and assembles them into a single global matrix by adding overlapping contributions at shared nodes.
        
        Parameters
        ----------
        local_matrix : np.ndarray
            Local element matrices with shape (n_elements, n_local_dofs, n_local_dofs)
            where n_local_dofs = polynomial_order + 1
            
        Returns
        -------
        np.ndarray
            Assembled global matrix with shape (n_global_dofs, n_global_dofs)
            where n_global_dofs = n_elements * polynomial_order + 1
            
        Notes
        -----
        The assembly process accounts for shared endpoints between adjacent elements,
        accumulating contributions from all elements that share a degree of freedom.
        """
        # Initialize global matrix
        global_size = self.finite_element_number * (self.physical_node_number - 1) + 1
        global_matrix = np.zeros((global_size, global_size), dtype=self.assembly_dtype)
        
        # Get assembly indices (mapping from local to global DOFs)
        rows_global, cols_global = self._build_assembly_indices()
        
        # Assemble: add local contributions to global matrix
        np.add.at(
            global_matrix,
            (rows_global, cols_global),
            np.asarray(local_matrix, dtype=self.assembly_dtype).reshape(-1),
        )
        
        return global_matrix
    
    
    def _build_assembly_indices_vector(self):
        """
        Return local→global indices for assembling element vectors into
        a global vector (with shared endpoints).
        
        This is similar to _build_assembly_indices but for vectors (1D).

        Primary FE space:
            local dofs per element       : N_grid = physical_node_number
            global stride (overlap by 1) : stride = N_grid - 1
            global dofs                  : N_elem * stride + 1
        
        Returns
        -------
        indices_global : np.ndarray
            1D array of global indices, shape (n_elements * n_local_dofs,)
            Maps local_vector.reshape(-1) to positions in global_vector
        
        Example
        -------
        For 3 elements with 3 nodes each (polynomial order 2):
            Element 0: local indices [0, 1, 2] → global indices [0, 1, 2]
            Element 1: local indices [0, 1, 2] → global indices [2, 3, 4]
            Element 2: local indices [0, 1, 2] → global indices [4, 5, 6]
        
        Returns: [0, 1, 2, 2, 3, 4, 4, 5, 6]
                  ^^^^^^^  ^^^^^^^  ^^^^^^^
                  elem 0   elem 1   elem 2
        
        Note the repeated indices (2, 2) and (4, 4) at element boundaries.
        """
        N_elem = self.finite_element_number  # e for element
        N_grid = self.physical_node_number       # g for grid

        # compute the assembly indices
        stride = N_grid - 1
        elem_array = np.arange(N_elem)           # [0, 1, 2, ...]
        grid_array = np.arange(N_grid)           # [0, 1, 2, ..., N_grid-1]
        
        # Global index for each (element, local_node) pair
        # I_eg[elem, local_idx] = global_idx
        I_eg = elem_array[:, None] * stride + grid_array[None, :]   # (N_elem, N_grid)
        
        # Flatten to 1D array
        indices_global = I_eg.reshape(-1)
        
        return indices_global
    
    
    def _build_assembly_indices(self):
        """
        Return clean local→global row/col indices for assembling element blocks into
        a global matrix (with shared endpoints).

        Primary FE space:
            local dofs per element       : N_g    = physical_node_number - 1
            global stride (overlap by 1) : stride = N_g - 1
            global dofs                  : N_e * stride + 1
        """
        # ----- primary FE space -----
        N_elem = self.number_of_finite_elements # e for element
        N_grid = self.physical_node_number      # g for grid

        # compute the assembly indices
        stride = N_grid - 1
        elem_array = np.arange(N_elem)
        grid_array = np.arange(N_grid)
        I_eg = elem_array[:, None]*stride + grid_array[None, :]   # (N_elem, N_grid)

        rows_primary = np.repeat(I_eg[:, :, None], N_grid, axis=2).reshape(-1)
        cols_primary = np.repeat(I_eg[:, None, :], N_grid, axis=1).reshape(-1)

        return rows_primary, cols_primary


    def _build_basis_on_arbitrary_grid_dict(self, given_grid: np.ndarray) -> Dict[int, np.ndarray]:
        """
        Evaluate Lagrange basis functions on arbitrary grid points, per element.

        Only elements that contain at least one point in ``given_grid`` appear in the
        returned dict (coarse uniform meshes may leave many elements out).

        Returns
        -------
        dict[int, np.ndarray]
            Keys are finite-element indices. Values have shape (1, n_points_in_element, n_basis).
        """

        # Internal helper: callers validate ``given_grid`` before calling.
        basis_on_arbitrary_grid_dict: Dict[int, np.ndarray] = {}
        for elem_idx in range(self.finite_element_number):
            idx_uniform = (
                (given_grid >= self.physical_nodes_reshaped[elem_idx:elem_idx + 1, 0])
                & (given_grid <= self.physical_nodes_reshaped[elem_idx:elem_idx + 1, -1])
            )
            if not np.any(idx_uniform):
                continue
            basis_func_in_current_element, _ = LagrangeShapeFunctions.lagrange_basis_and_derivatives(
                x_node = self.physical_nodes_reshaped[elem_idx:elem_idx + 1, :],
                x_eval = given_grid[idx_uniform],
            )
            basis_on_arbitrary_grid_dict[elem_idx] = basis_func_in_current_element
        return basis_on_arbitrary_grid_dict


    @property
    def H_kinetic(self) -> np.ndarray:
        return self.get_H_kinetic()

    @property
    def H_r_inv_sq(self) -> np.ndarray:
        return self.get_H_r_inv_sq()

    @property
    def S(self) -> np.ndarray:
        return self.get_S()

    @property
    def S_inv_sqrt(self) -> np.ndarray:
        return self.get_S_inv_sqrt()

    @property
    def laplacian(self) -> np.ndarray:
        return self.get_laplacian()

    @property
    def derivative_matrix(self) -> np.ndarray:
        # return self.get_derivative_matrix()
        # Now using the quadrature basis as the basis nodes
        return self.get_derivative_matrix_with_quadrature_basis()
        
    @property
    def derivative_matrix_with_quadrature_basis(self) -> np.ndarray:
        return self.get_derivative_matrix_with_quadrature_basis()

    @property
    def global_interpolation_matrix(self) -> np.ndarray:
        """
        Property accessor for global interpolation matrix.
        Used for HF/PBE0 calculations.
        """
        return self.get_global_interpolation_matrix()
    
    @property
    def lagrange_basis_pseudoinverse(self) -> np.ndarray:
        return self.get_lagrange_basis_pseudoinverse()



    @property
    def grid_data(self) -> GridData:
        """
        Property accessor for grid data.
        Used for response function calculation, etc.
        """
        return GridData(
            finite_element_number = self.finite_element_number,
            physical_nodes        = self.physical_nodes,
            quadrature_nodes      = self.quadrature_nodes,
            quadrature_weights    = self.quadrature_weights,
        )


    def evaluate_orbitals_on_arbitrary_grid(
        self,
        given_grid                   : np.ndarray,
        orbital_coefficients         : np.ndarray,
        basis_on_arbitrary_grid_dict : Optional[Dict[int, np.ndarray]] = None,
    ) -> np.ndarray:
        """
        Interpolate FE nodal orbital coefficients onto an arbitrary grid.

        Use when SCF eigenvectors are stored as global FE node coefficients
        (``symmetrize=False`` path). For quadrature-point orbitals, use
        ``evaluate_quantites_on_arbitrary_grid`` instead.

        Parameters
        ----------
        given_grid : np.ndarray, shape (n_points,)
            Monotonically increasing points within the physical domain.
        orbital_coefficients : np.ndarray
            Shape ``(n_global_fe_dofs - 1, n_orbitals)`` before internal padding.
        basis_on_arbitrary_grid_dict : dict[int, np.ndarray], optional
            Precomputed basis from ``_build_basis_on_arbitrary_grid_dict``.

        Returns
        -------
        np.ndarray, shape (n_points, n_orbitals)
            Orbitals sampled on ``given_grid``.
        """
        # Validate input: grid must be monotonically increasing
        assert np.all(np.diff(given_grid) > 0.0)

        # Validate input: grid must be within physical domain
        assert np.all(given_grid >= self.physical_nodes[0]) and \
               np.all(given_grid <= self.physical_nodes[-1])
            
        # Validate input: FE nodal coefficients (unpadded global dof count)
        assert orbital_coefficients.shape[0] == (self.finite_element_number * (self.physical_node_number - 1) - 1)

        if basis_on_arbitrary_grid_dict is None:
            basis_on_arbitrary_grid_dict = self._build_basis_on_arbitrary_grid_dict(given_grid)


        orbital_coefficients = np.pad(orbital_coefficients, ((1, 1), (0, 0)))
        n_elem = self.finite_element_number
        n_basis = self.physical_node_number

        orbitals_arbitrary_grid = np.zeros((len(given_grid), orbital_coefficients.shape[1]))
        for elem_idx in range(n_elem):
            idx_FE = np.arange(elem_idx * (n_basis - 1), (elem_idx + 1) * (n_basis - 1) + 1)
            idx_uniform = (
                (given_grid >= self.physical_nodes_reshaped[elem_idx:elem_idx + 1, 0])
                & (given_grid <= self.physical_nodes_reshaped[elem_idx:elem_idx + 1, -1])
            )
            if not np.any(idx_uniform):
                continue
            orbitals_arbitrary_grid[idx_uniform, :] = (
                basis_on_arbitrary_grid_dict[elem_idx][0, :, :] @ orbital_coefficients[idx_FE, :]
            )
        return orbitals_arbitrary_grid


    def evaluate_quantites_on_arbitrary_grid(
        self,
        given_grid                   : np.ndarray,
        field_values                 : np.ndarray,
        basis_on_arbitrary_grid_dict : Optional[Dict[int, np.ndarray]] = None,
    ) -> np.ndarray:
        """
        Evaluate quadrature-point fields on an arbitrary grid.

        For each finite element, quadrature values are mapped to nodal coefficients
        via the cached Lagrange-basis pseudoinverse, then interpolated to
        ``given_grid`` using precomputed basis values on that grid.

        Parameters
        ----------
        given_grid : np.ndarray, shape (n_points,)
            Monotonically increasing evaluation points within the physical domain.
        field_values : np.ndarray, shape (n_elem * n_quad, n_fields)
            Field values at global quadrature points. Use ``field[:, None]`` for a
            single field.
        basis_on_arbitrary_grid_dict : dict[int, np.ndarray], optional
            Precomputed Lagrange basis on ``given_grid``, keyed by element index.
            Each value has shape ``(1, n_points_in_element, n_basis)``.
            If ``None``, the basis is built internally for this call.
            Pass the dict from ``_build_basis_on_arbitrary_grid_dict`` to avoid
            rebuilding when evaluating multiple fields on the same grid.

        Returns
        -------
        np.ndarray, shape (n_points, n_fields)
            Field values interpolated onto ``given_grid``.
        """
        # Validate input: grid must be monotonically increasing
        assert np.all(np.diff(given_grid) > 0.0)

        # Validate input: grid must be within physical domain
        assert np.all(given_grid >= self.physical_nodes[0]) and \
               np.all(given_grid <= self.physical_nodes[-1])
        
        # Validate input: field values must be a 1D or 2D array with shape (n_elem * n_quad, n_fields)
        if field_values.ndim == 1:
            field_values = field_values[:, None]
        assert field_values.ndim == 2
        assert field_values.shape[0] == self.finite_element_number * self.quadrature_node_number

        n_elem = self.finite_element_number
        n_quad = self.quadrature_node_number
        basis_pseudoinverse = self.lagrange_basis_pseudoinverse

        if basis_on_arbitrary_grid_dict is None:
            basis_on_arbitrary_grid_dict = self._build_basis_on_arbitrary_grid_dict(given_grid)

        quantities_arbitrary_grid = np.zeros((len(given_grid), field_values.shape[1]))
        for elem_idx in range(n_elem):
            idx_uniform = (
                (given_grid >= self.physical_nodes_reshaped[elem_idx:elem_idx + 1, 0])
                & (given_grid <= self.physical_nodes_reshaped[elem_idx:elem_idx + 1, -1])
            )
            if not np.any(idx_uniform):
                continue
            field_values_coeffs_from_pseudoinverse = (
                basis_pseudoinverse[elem_idx, :, :]
                @ field_values[elem_idx * n_quad:(elem_idx + 1) * n_quad, :]
            )
            quantities_arbitrary_grid[idx_uniform, :] = (
                basis_on_arbitrary_grid_dict[elem_idx] @ field_values_coeffs_from_pseudoinverse
            )

        return quantities_arbitrary_grid


    def evaluate_single_field_on_grid(
        self,
        given_grid                   : np.ndarray,
        field_values                 : np.ndarray,
        basis_on_arbitrary_grid_dict : Optional[Dict[int, np.ndarray]] = None,
    ) -> np.ndarray:
        """
        Evaluate a single field on a given grid using Lagrange interpolation.
        
        This function takes a field represented by its values at quadrature
        points and evaluates it at arbitrary grid points using finite element
        Lagrange basis functions. For each grid point, the function:
        1. Identifies which finite element contains the point
        2. Converts quadrature point values to nodal coefficients (using pseudoinverse)
        3. Evaluates the Lagrange basis functions at that point
        4. Computes the orbital value as a linear combination of basis functions
        
        Parameters
        ----------
        given_grid : np.ndarray, shape (n_points,)
            Grid points where the field should be evaluated.
            Must be monotonically increasing and within the physical domain.
        field_values : np.ndarray, shape (n_elem * n_quad,)
            Field values at global quadrature points.
        basis_on_arbitrary_grid_dict : dict[int, np.ndarray], optional
            Precomputed basis from ``_build_basis_on_arbitrary_grid_dict``.
            Reuse when mapping several fields onto the same ``given_grid``.

        Returns
        -------
        field_on_grid : np.ndarray, shape (n_points,)
            Field values on ``given_grid``.

        Notes
        -----
        Delegates to ``evaluate_quantites_on_arbitrary_grid`` (vectorized per
        element). Input is quadrature-point data, not FE nodal coefficients.
        
        Example
        -------
        >>> # Given field values at quadrature points
        >>> field_values = np.array([...])  # shape: (n_elem * n_quad,)
        >>> # Evaluate on a uniform grid
        >>> uniform_grid = np.linspace(0, domain_size, 1000)
        >>> field_on_grid = ops_builder.evaluate_single_field_on_grid(
        ...     given_grid=uniform_grid,
        ...     field_values=field_values
        ... )
        >>> # field_on_grid.shape = (1000,)
        """
        # Validate input: field values must be a 1D array with shape (n_elem * n_quad,)
        assert field_values.ndim == 1
        assert field_values.shape[0] == self.finite_element_number * self.quadrature_node_number
        
        return self.evaluate_quantites_on_arbitrary_grid(
            given_grid                   = given_grid,
            field_values                 = field_values,
            basis_on_arbitrary_grid_dict = basis_on_arbitrary_grid_dict,
        ).reshape(-1)



    @staticmethod
    def build_cubic_spline_derivative_matrix(
        r                : np.ndarray,
        left_derivative  : float = 0.0,
        right_derivative : float = 0.0,
    ) -> np.ndarray:
        """
        Build a 2D derivative matrix for a 1D cubic spline on ``r``.

        For grid values ``f`` on ``r``, ``D @ f`` gives ``df/dr`` at each grid point,
        where ``D[i, j] = ∂f/∂f_j`` evaluated at ``r_i``. The spline uses clamped
        boundary conditions ``f'(r_0) = left_derivative`` and
        ``f'(r_{n-1}) = right_derivative`` (default: zero slope at both ends).

        Parameters
        ----------
        r : np.ndarray
            Strictly increasing 1D radial grid, shape ``(n,)``.
        left_derivative, right_derivative : float
            First-derivative boundary values at the left and right endpoints.

        Returns
        -------
        np.ndarray
            Shape ``(n, n)``, ``float64``.
        """
        from scipy.interpolate import CubicSpline

        r = np.asarray(r, dtype=np.float64).reshape(-1)
        n = r.size
        if n < 2:
            raise RuntimeError(
                "build_cubic_spline_derivative_matrix requires at least 2 grid points, "
                f"got {n}."
            )
        if not np.all(np.diff(r) > 0.0):
            raise RuntimeError("Grid r must be strictly increasing.")

        bc_type = ((1, float(left_derivative)), (1, float(right_derivative)))
        D = np.zeros((n, n), dtype=np.float64)
        eye = np.eye(n, dtype=np.float64)
        for j in range(n):
            spline = CubicSpline(r, eye[:, j], bc_type=bc_type)
            D[:, j] = spline(r, nu=1)
        return D





def f(r: np.ndarray) -> np.ndarray:
    """Benchmark function: f(r) = sin(r) + 0.35 * sin(2r)."""
    return np.sin(r) + 0.35 * np.sin(2.0 * r)


def f_prime(r: np.ndarray) -> np.ndarray:
    """Analytical derivative: f'(r) = cos(r) + 0.70 * cos(2r)."""
    return np.cos(r) + 0.70 * np.cos(2.0 * r)


"""Compute numerical/analytical derivatives for a smooth periodic benchmark."""
grid_data = GridData.from_basic(
    domain_size = 4.0 * np.pi,
    finite_element_number = 6,
    polynomial_order = 20,
    quadrature_point_number = 24,
    mesh_type = "exponential",
    mesh_concentration = 2.0, 
)

ops_builder = RadialOperatorsBuilder.from_grid_data(grid_data, verbose=False)

D = ops_builder.get_derivative_matrix_with_quadrature_basis()
r = ops_builder.quadrature_nodes_reshaped

f_values = f(r)
f_ref_prime = f_prime(r)

f_num_prime = np.einsum("eij,ej->ei", D, f_values)
abs_error = np.abs(f_num_prime - f_ref_prime)
max_abs_error = np.max(abs_error)
