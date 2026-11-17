"""Filter and projection: the properties docs/topopt.md relies on."""

import numpy as np

from adjoint_helmholtz.topopt import binary_fraction, conic_filter, project, project_derivative


def test_filter_averages():
    # Rows sum to 1 with non-negative weights, so a uniform density stays uniform and 0 <= W rho <= 1.
    W = conic_filter((12, 9), radius_px=2.5)
    assert np.allclose(W.sum(axis=1), 1)
    assert W.min() >= 0
    assert np.allclose(W @ np.full(12 * 9, 0.3), 0.3)


def test_filter_reaches_radius_and_no_further():
    nx, ny, radius = 15, 15, 3.0
    W = conic_filter((nx, ny), radius)
    centre = 7 * ny + 7
    row = W[[centre]].toarray().reshape(nx, ny)
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    distance = np.hypot(i - 7, j - 7)
    assert np.all(row[distance >= radius] == 0)
    assert np.all(row[distance < radius] > 0)


def test_projection_end_points_and_threshold():
    for beta in [1, 8, 64]:
        assert np.isclose(project(0.0, beta), 0)
        assert np.isclose(project(1.0, beta), 1)
        assert np.isclose(project(0.5, beta), 0.5)
    # Large beta approaches a step at eta = 0.5.
    assert project(0.45, 128) < 0.01 and project(0.55, 128) > 0.99


def test_projection_derivative_matches_central_differences():
    rho_t = np.linspace(0, 1, 23)
    step = 1e-6
    for beta in [1, 8, 64]:
        fd = (project(rho_t + step, beta) - project(rho_t - step, beta)) / (2 * step)
        assert np.allclose(project_derivative(rho_t, beta), fd, rtol=1e-6, atol=1e-8)


def test_binary_fraction():
    assert binary_fraction(np.array([0.0, 0.04, 0.5, 0.96, 1.0])) == 0.8