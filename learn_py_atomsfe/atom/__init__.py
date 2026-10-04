"""
SPARC-atomSFE Atomic DFT Solver Package

This package provides a comprehensive implementation of Atomic Density Functional Theory
solver using finite element method.

Main class:
    AtomicDFTSolver: Main solver class for atomic DFT calculations

Example:
    >>> from atom import AtomicDFTSolver
    >>> solver = AtomicDFTSolver(atomic_number=13, xc_functional="GGA_PBE")
    >>> results = solver.solve()
"""

__version__ = "0.1.0"


import sys
import os

# ---------------------------------------------------------------------------
# Record whether NumPy was already imported BEFORE this package was imported.
# This is the key information for determining BLAS thread safety later.
# ---------------------------------------------------------------------------
_NUMPY_IMPORTED_BEFORE_ATOMIC = ("numpy" in sys.modules)


# ---------------------------------------------------------------------------
# Record whether BLAS environment variables already indicate single-threaded
# execution. This allows us to avoid warnings even when NumPy was imported
# earlier, because the user may have correctly configured BLAS beforehand.
# ---------------------------------------------------------------------------

def _env_says_single_thread() -> bool:
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        val = os.environ.get(var)
        if val is None:
            return False
        if val.strip() not in ("1", "1.0"):
            return False
    return True

_BLAS_ENV_SINGLE_THREADED = _env_says_single_thread()

# ---------------------------------------------------------------------------
# Record whether threadpoolctl is installed.
# ---------------------------------------------------------------------------
def _threadpoolctl_installed() -> bool:
    try:
        import threadpoolctl
        return True
    except ImportError:
        return False

_THREADPOOLCTL_INSTALLED = _threadpoolctl_installed()

# BLAS/OpenMP thread counts are NOT forced here: whatever the environment says is
# what the run uses.  _BLAS_ENV_SINGLE_THREADED and _THREADPOOLCTL_INSTALLED are
# still recorded above, for the parallel sections to report on.



# Lazy import: solver is only loaded when AtomicDFTSolver is actually used.
# (ML-XC is not wired through the public solver API in this fork, so importing the solver
# no longer pulls torch via ml_xc unless lower-level modules import it.)
__all__ = ["AtomicDFTSolver"]


def __getattr__(name: str):
    if name == "AtomicDFTSolver":
        from .solver import AtomicDFTSolver
        return AtomicDFTSolver
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
