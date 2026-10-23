Unit tests for derivative-matrix accuracy in radial mesh operators.

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
