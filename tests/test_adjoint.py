"""Adjoint gradient against central finite differences (docs/adjoint.md)."""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from adjoint_helmholtz.adjoint import intensity_and_gradient
from adjoint_helmholtz.helmholtz1d import helmholtz_matrix, point_source

K0, H, N_PML, N = 2 * np.pi, 1 / 20, 10, 101


def setup():
    """A random slab between a point source (node 20) and a probe (node 80)."""
    rng = np.random.default_rng(0)
    eps_r = np.ones(N)
    eps_r[40:60] = 2 + 0.3 * rng.random(20)
    b = point_source(N, 20, H)
    c = np.zeros(N)
    c[80] = 1
    return eps_r, b, c


def intensity(eps_r, b, c):
    u = spsolve(helmholtz_matrix(eps_r, H, K0, N_PML), b)
    return abs(np.vdot(c, u)) ** 2


def test_gradient_matches_central_differences():
    eps_r, b, c = setup()
    _, grad = intensity_and_gradient(helmholtz_matrix(eps_r, H, K0, N_PML), b, c, K0)
    step = 1e-4
    # Nodes 5 and 95 are in the PML: a wrong transpose (A instead of A^T) only shows there.
    for i in [5, 30, 45, 55, 70, 95]:
        plus, minus = eps_r.copy(), eps_r.copy()
        plus[i] += step
        minus[i] -= step
        fd = (intensity(plus, b, c) - intensity(minus, b, c)) / (2 * step)
        assert abs(fd - grad[i]) < 1e-6 * np.max(np.abs(grad)), i


def test_adjoint_equals_direct_sensitivity_3x3():
    # Smallest case: a random non-symmetric complex 3x3 A. The direct way needs one solve per
    # parameter, du_i = -A^{-1} (dA/d eps_i) u, and uses no transpose. Both must agree to rounding error.
    rng = np.random.default_rng(1)
    k0 = 1.3
    A0 = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    eps_r = 1 + rng.random(3)
    A = A0 + k0**2 * np.diag(eps_r)
    b = rng.normal(size=3) + 1j * rng.normal(size=3)
    c = rng.normal(size=3) + 1j * rng.normal(size=3)

    J, grad = intensity_and_gradient(sp.csc_matrix(A), b, c, k0)

    u = np.linalg.solve(A, b)
    a = np.vdot(c, u)
    direct = []
    for i in range(3):
        dA = np.zeros((3, 3))
        dA[i, i] = k0**2  # docs/adjoint.md, Eq. 1
        du = -np.linalg.solve(A, dA @ u)
        direct.append(2 * np.real(np.conj(a) * np.vdot(c, du)))  # dJ = 2 Re(conj(a) da), da = c^H du
    assert np.isclose(J, abs(a) ** 2)
    assert np.allclose(grad, direct, rtol=1e-12, atol=0)