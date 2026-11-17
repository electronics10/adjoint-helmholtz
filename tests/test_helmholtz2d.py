"""Checks of the 2D solver: Hankel function, convergence, PML reflection, symmetry of S A."""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.special import hankel1

from adjoint_helmholtz.helmholtz2d import helmholtz_matrix, point_source
from adjoint_helmholtz.pml import stretch_factors

K0 = 2 * np.pi  # wavelength 1 (any length unit)


def solve_point_source(h, n_phys, n_pml):
    """Vacuum, square grid, unit point source in the middle. Returns u (n, n), r (n, n) and the physical mask."""
    n = n_phys + 2 * n_pml
    c = n // 2
    A = helmholtz_matrix(np.ones((n, n)), h, K0, n_pml)
    u = spsolve(A, point_source((n, n), c, c, h)).reshape(n, n)
    x = (np.arange(n) - c) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    physical = np.zeros((n, n), bool)
    physical[n_pml : n - n_pml, n_pml : n - n_pml] = True
    return u, np.hypot(X, Y), physical


def relative_error(u, reference):
    return np.max(np.abs(u - reference)) / np.max(np.abs(reference))


def hankel_error(h):
    """Error against G = (i/4) H0(k r) in the physical box, at least half a wavelength from the source."""
    u, r, physical = solve_point_source(h, n_phys=round(2 / h) + 1, n_pml=round(0.5 / h))
    far = physical & (r >= 0.5)
    G = 0.25j * hankel1(0, K0 * r[far])
    return relative_error(u[far], G)


def test_matches_hankel_function():
    assert hankel_error(1 / 40) < 0.01


def test_second_order_convergence():
    errors = [hankel_error(h) for h in [1 / 20, 1 / 40, 1 / 80]]
    orders = np.log2(np.array(errors[:-1]) / np.array(errors[1:]))
    assert np.all((orders > 1.9) & (orders < 2.1)), orders


def pml_error(n_pml, h=1 / 20, n_phys=41, n_pad=80):
    """Error in the physical box against the same box cut out of a much larger domain (PML far away)."""
    u, _, physical = solve_point_source(h, n_phys, n_pml)
    big, _, _ = solve_point_source(h, n_phys + 2 * n_pad, n_pml=20)
    inner = slice(20 + n_pad, 20 + n_pad + n_phys)
    return relative_error(u[physical].reshape(n_phys, n_phys), big[inner, inner])


def test_pml_reflection_drops_with_thickness():
    e10, e20 = pml_error(10), pml_error(20)
    assert e10 < 1e-4
    assert e20 < e10 / 5


def test_scaled_matrix_is_symmetric():
    # S = diag(s_x(i) s_y(j)) makes A complex symmetric (docs/fdfd2d.md, Eq. 3).
    rng = np.random.default_rng(0)
    nx, ny, h, n_pml = 30, 24, 1 / 20, 8
    eps_r = 1 + rng.random((nx, ny)) + 0.1j * rng.random((nx, ny))
    A = helmholtz_matrix(eps_r, h, K0, n_pml)
    sx, _ = stretch_factors(nx, n_pml, h, K0)
    sy, _ = stretch_factors(ny, n_pml, h, K0)
    SA = sp.diags(np.outer(sx, sy).ravel()) @ A
    assert abs(SA - SA.T).max() < 1e-12 * abs(SA).max()