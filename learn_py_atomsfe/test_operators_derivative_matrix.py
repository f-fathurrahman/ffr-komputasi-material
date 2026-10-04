"""Unit tests for derivative-matrix accuracy in radial mesh operators.

This module validates that the derivative matrix constructed by
`RadialOperatorsBuilder` correctly approximates the first derivative of a
known smooth function at element-wise quadrature points.

Benchmark function used throughout this file:
    f(r) = sin(r) + 0.35 * sin(2r)
    f'(r) = cos(r) + 0.70 * cos(2r)

What is tested:
1) Build a radial grid on domain [0, 4*pi] using `GridData.from_basic`.
2) Construct operators with `RadialOperatorsBuilder.from_grid_data`.
3) Compute derivative matrix D with quadrature-point basis.
4) Evaluate f(r) on quadrature nodes and compute numerical derivative D @ f.
5) Compare numerical derivative against analytical f'(r) and assert max-abs
   error is below tolerance.

Script mode (`python tests/unit/test_operators_derivative_matrix.py`) also:
- saves a full-domain comparison plot with a secondary y-axis showing
  absolute error.

How to run:
- Run only this file with pytest:
  python -m pytest tests/unit/test_operators_derivative_matrix.py -q

- Run script mode for visualization:
  python tests/unit/test_operators_derivative_matrix.py
"""

from pathlib import Path

import numpy as np

from atom.mesh.operators import GridData, RadialOperatorsBuilder


def f(r: np.ndarray) -> np.ndarray:
    """Benchmark function: f(r) = sin(r) + 0.35 * sin(2r)."""
    return np.sin(r) + 0.35 * np.sin(2.0 * r)


def f_prime(r: np.ndarray) -> np.ndarray:
    """Analytical derivative: f'(r) = cos(r) + 0.70 * cos(2r)."""
    return np.cos(r) + 0.70 * np.cos(2.0 * r)


def _compute_derivative_matrix_benchmark():
    """Compute numerical/analytical derivatives for a smooth periodic benchmark."""
    grid_data = GridData.from_basic(
        domain_size             = 4.0 * np.pi,
        finite_element_number   = 6,
        polynomial_order        = 20,
        quadrature_point_number = 24,
        mesh_type               = "exponential",
        mesh_concentration      = 2.0, 
    )
    ops_builder = RadialOperatorsBuilder.from_grid_data(grid_data, verbose=False)
    D = ops_builder.get_derivative_matrix_with_quadrature_basis()
    r = ops_builder.quadrature_nodes_reshaped
    f_values = f(r)
    f_ref_prime = f_prime(r)
    f_num_prime = np.einsum("eij,ej->ei", D, f_values)
    abs_error = np.abs(f_num_prime - f_ref_prime)
    max_abs_error = np.max(abs_error)
    return r, f_ref_prime, f_num_prime, abs_error, max_abs_error


def test_derivative_matrix_against_analytic_periodic_derivative():
    """Compare Df against analytic f'(r) on quadrature points.

    The test builds a radial operators object, computes the derivative matrix,
    applies it to smooth periodic values, and checks the numerical derivative
    against the exact analytical derivative.
    """
    _, _, _, _, max_abs_error = _compute_derivative_matrix_benchmark()
    assert max_abs_error < 1e-3


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    r, f_ref_prime, f_num_prime, abs_error, max_abs_error = _compute_derivative_matrix_benchmark()
    print(f"max_abs_error: {max_abs_error:.6e}")

    r_flat = r.reshape(-1)
    ref_flat = f_ref_prime.reshape(-1)
    num_flat = f_num_prime.reshape(-1)
    err_flat = abs_error.reshape(-1)
    sort_idx = np.argsort(r_flat)

    fig, ax1 = plt.subplots(figsize=(9, 5))
    ax1.plot(r_flat[sort_idx], ref_flat[sort_idx], label="Analytical df/dr", linewidth=2.0, color="C0")
    ax1.plot(r_flat[sort_idx], num_flat[sort_idx], "--", label="Numerical D@f", linewidth=1.6, color="C1")
    ax1.set_xlabel("radius (Bohr)")
    ax1.set_ylabel("Derivative value (analytical and numerical), $df/dr$")
    ax1.set_title("Derivative matrix vs analytical derivative (domain: 0 to 4π)")
    ax1.grid(True, alpha=0.3)

    formula_text = (
        r"$f(r)=\sin(r)+0.35\sin(2r)$" "\n"
        r"$f'(r)=\cos(r)+0.70\cos(2r)$"
    )
    ax1.text(
        0.18,
        0.90,
        formula_text,
        transform=ax1.transAxes,
        va="top",
        ha="left",
        fontsize=13,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="0.6"),
    )

    ax2 = ax1.twinx()
    ax2.plot(r_flat[sort_idx], err_flat[sort_idx], label="Absolute error", linewidth=1.4, color="C3", alpha=0.45)
    ax2.set_ylabel("absolute error", color="C3")
    ax2.tick_params(axis="y", labelcolor="C3")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    fig.tight_layout()
    output_path = Path(__file__).resolve().parent / "check_derivative_matrix.png"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(output_path, dpi=200)
    plt.show()
    plt.close(fig)
