"""Stretched-coordinate PML profile along one axis (docs/fdfd1d.md, Eqs. 3-4)."""

import numpy as np


def stretch_factors(n, n_pml, h, k, reflection=1e-8, order=3):
    """Stretch factors s = 1 + i*alpha on a line of n nodes x_i = i*h.

    Returns s_node, shape (n,), at x_i, and s_half, shape (n + 1,), at x_{i - 1/2}.
    The PML fills the outermost n_pml cells on each side; k is the wavenumber it is tuned for.
    """
    s_node = np.ones(n, complex)
    s_half = np.ones(n + 1, complex)
    if n_pml == 0:
        return s_node, s_half

    thickness = n_pml * h
    alpha_max = (order + 1) * np.log(1 / reflection) / (2 * k * thickness)  # Eq. 4
    x_left = n_pml * h  # where the PML starts, seen from the inside
    x_right = (n - 1 - n_pml) * h

    def alpha(x):
        depth = np.maximum(x_left - x, 0) + np.maximum(x - x_right, 0)
        depth = np.minimum(depth, thickness)
        return alpha_max * (depth / thickness) ** order  # Eq. 3

    x_node = np.arange(n) * h
    x_half = (np.arange(n + 1) - 0.5) * h
    return 1 + 1j * alpha(x_node), 1 + 1j * alpha(x_half)