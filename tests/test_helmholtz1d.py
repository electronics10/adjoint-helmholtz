"""Checks of the 1D solver that need no trusted code: known answers, convergence, symmetry."""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from adjoint_helmholtz.helmholtz1d import helmholtz_matrix, point_source
from adjoint_helmholtz.pml import stretch_factors

K0 = 2 * np.pi  # wavelength 1 (any length unit)


def solve_point_source(h, n_phys, n_pml):
    """Vacuum, unit point source in the middle. Returns x, u, the source index and the physical slice."""
    n = n_phys + 2 * n_pml
    j = n // 2
    A = helmholtz_matrix(np.ones(n), h, K0, n_pml)
    u = spsolve(A, point_source(n, j, h))
    x = np.arange(n) * h
    return x, u, j, slice(n_pml, n - n_pml)


def relative_error(u, reference):
    return np.max(np.abs(u - reference)) / np.max(np.abs(reference))


def test_matches_continuous_green_function():
    # G = i/(2k) exp(ik|x - x0|), the response to -delta (conventions.md, item 2)
    x, u, j, phys = solve_point_source(h=1 / 40, n_phys=81, n_pml=20)
    G = 1j / (2 * K0) * np.exp(1j * K0 * np.abs(x - x[j]))
    assert relative_error(u[phys], G[phys]) < 0.01


def test_second_order_convergence():
    # Fixed physical box (2 wavelengths) and PML thickness (half a wavelength); halve h three times.
    errors = []
    for h in [1 / 20, 1 / 40, 1 / 80, 1 / 160]:
        x, u, j, phys = solve_point_source(h, n_phys=round(2 / h) + 1, n_pml=round(0.5 / h))
        G = 1j / (2 * K0) * np.exp(1j * K0 * np.abs(x - x[j]))
        errors.append(relative_error(u[phys], G[phys]))
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all((orders > 1.9) & (orders < 2.1)), orders


def pml_error(n_pml, h=1 / 20):
    """Error against the exact solution of the infinite grid, so only the PML is measured.

    Infinite grid: u_i = -h q^|i-j| / (2i sin Kh), q = exp(iKh), cos Kh = 1 - (k0 h)^2 / 2
    (docs/fdfd1d.md, Eq. 6).
    """
    x, u, j, phys = solve_point_source(h, n_phys=81, n_pml=n_pml)
    Kh = np.arccos(1 - (K0 * h) ** 2 / 2)
    i = np.arange(len(u))
    exact = -h * np.exp(1j * Kh * np.abs(i - j)) / (2j * np.sin(Kh))
    return relative_error(u[phys], exact[phys])


def test_pml_reflection_drops_with_thickness():
    e10, e20 = pml_error(10), pml_error(20)
    assert e10 < 1e-3
    assert e20 < e10 / 5


def test_scaled_matrix_is_symmetric():
    # Reciprocity: S A is complex symmetric, S = diag(s_node), even with loss (conventions.md, item 5).
    rng = np.random.default_rng(0)
    n, h, n_pml = 60, 1 / 20, 10
    eps_r = 1 + rng.random(n) + 0.1j * rng.random(n)
    A = helmholtz_matrix(eps_r, h, K0, n_pml)
    s, _ = stretch_factors(n, n_pml, h, K0)
    SA = sp.diags(s) @ A
    assert abs(SA - SA.T).max() < 1e-12 * abs(SA).max()
    assert abs(A - A.T).max() > 1e-3 * abs(A).max()  # A itself is not symmetric
