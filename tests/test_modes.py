"""Mode solver, one-way mode source and monitor, in a straight silicon waveguide."""

import numpy as np
from scipy.sparse.linalg import spsolve

from adjoint_helmholtz.helmholtz2d import helmholtz_matrix
from adjoint_helmholtz.modes import discrete_beta, solve_modes
from adjoint_helmholtz.ports import incident_field, mode_monitor, mode_source, total_field_mask

WAVELENGTH_UM = 1.55
K0 = 2 * np.pi / WAVELENGTH_UM
H_UM = 0.025
N_PML = 40  # 1 um


def straight_guide(width_um=0.4, length_um=4.0, height_um=3.0):
    """Si core (eps 12.1) in SiO2 (eps 2.1), along x. Returns eps_r (nx, ny)."""
    nx, ny = round(length_um / H_UM) + 1, round(height_um / H_UM) + 1
    y = np.arange(ny) * H_UM
    eps_line = np.where(np.abs(y - y[-1] / 2) <= width_um / 2, 12.1, 2.1)
    return np.tile(eps_line, (nx, 1))


def test_mode_is_guided():
    eps_r = straight_guide()
    beta2, _ = solve_modes(eps_r[0], H_UM, K0, N_PML, num_modes=2)
    n_eff = np.sqrt(beta2.real) / K0
    assert np.all((n_eff > np.sqrt(2.1)) & (n_eff < np.sqrt(12.1)))
    assert n_eff[0] > n_eff[1]


def test_source_is_one_way_and_transmission_is_one():
    eps_r = straight_guide()
    shape = eps_r.shape
    A = helmholtz_matrix(eps_r, H_UM, K0, N_PML)
    beta2, e = solve_modes(eps_r[0], H_UM, K0, N_PML)
    beta_d = discrete_beta(beta2[0], H_UM)
    i_src, i_mon = N_PML + 10, shape[0] - N_PML - 10

    u = spsolve(A, mode_source(A, shape, H_UM, e[:, 0], beta_d, axis=0, index=i_src))

    # Outside the PML, the field is exactly the launched mode on one side and zero on the other.
    expected = total_field_mask(shape, 0, i_src, +1) * incident_field(shape, H_UM, e[:, 0], beta_d, 0, i_src, +1)
    physical = (slice(N_PML, shape[0] - N_PML), slice(N_PML, shape[1] - N_PML))
    error = u.reshape(shape) - expected
    assert np.max(np.abs(error[physical])) < 1e-5 * np.max(np.abs(expected))

    c = mode_monitor(shape, e[:, 0], axis=0, index=i_mon)
    assert abs(abs(np.vdot(c, u)) ** 2 - 1) < 1e-5