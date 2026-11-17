"""2D adjoint gradient of a mode transmission against central finite differences (docs/adjoint.md)."""

import numpy as np
from scipy.sparse.linalg import spsolve

from adjoint_helmholtz.adjoint import intensity_and_gradient
from adjoint_helmholtz.helmholtz2d import helmholtz_matrix
from adjoint_helmholtz.modes import discrete_beta, solve_modes
from adjoint_helmholtz.ports import mode_monitor, mode_source

K0 = 2 * np.pi / 1.55
H_UM = 0.04
N_PML = 15


def scattering_problem():
    """A Si guide along x, cut by a random 'design' block; mode source on the left, monitor on the right.

    Returns eps_r (nx, ny), b, c and the design-block slices.
    """
    nx, ny = round(3.0 / H_UM) + 1 + 2 * N_PML, round(2.4 / H_UM) + 1 + 2 * N_PML
    y = np.arange(ny) * H_UM
    eps_r = np.tile(np.where(np.abs(y - y[-1] / 2) <= 0.2, 12.1, 2.1), (nx, 1))
    design = (slice(N_PML + 20, nx - N_PML - 20), slice(ny // 2 - 15, ny // 2 + 15))
    rng = np.random.default_rng(1)
    eps_r[design] = 2.1 + 10 * rng.random(eps_r[design].shape)

    # The mode comes from the undisturbed input guide (column 0). b and c stay fixed when eps_r changes.
    beta2, e = solve_modes(eps_r[0], H_UM, K0, N_PML)
    A = helmholtz_matrix(eps_r, H_UM, K0, N_PML)
    b = mode_source(A, eps_r.shape, H_UM, e[:, 0], discrete_beta(beta2[0], H_UM), axis=0, index=N_PML + 5)
    c = mode_monitor(eps_r.shape, e[:, 0], axis=0, index=nx - N_PML - 5)
    return eps_r, b, c, design


def transmission(eps_r, b, c):
    u = spsolve(helmholtz_matrix(eps_r, H_UM, K0, N_PML), b)
    return abs(np.vdot(c, u)) ** 2


def test_gradient_matches_central_differences_on_20_pixels():
    eps_r, b, c, design = scattering_problem()
    nx, ny = eps_r.shape
    _, grad = intensity_and_gradient(helmholtz_matrix(eps_r, H_UM, K0, N_PML), b, c, K0)
    grad = grad.reshape(nx, ny)

    rng = np.random.default_rng(2)
    rows = rng.integers(design[0].start, design[0].stop, 16)
    cols = rng.integers(design[1].start, design[1].stop, 16)
    pixels = list(zip(rows, cols))
    # Four pixels in the PML, one per side: a wrong transpose (A instead of A^T) is
    # below 1e-7 inside the design block and only shows here (docs/adjoint.md).
    pixels += [(3, ny // 2), (nx - 4, ny // 2 + 2), (nx // 2, 3), (nx // 2, ny - 4)]

    step = 1e-3
    for p in pixels:
        plus, minus = eps_r.copy(), eps_r.copy()
        plus[p] += step
        minus[p] -= step
        fd = (transmission(plus, b, c) - transmission(minus, b, c)) / (2 * step)
        assert abs(fd - grad[p]) < 1e-6 * np.max(np.abs(grad)), p