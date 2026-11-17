"""Sparse finite-difference operators with stretched-coordinate PML (docs/fdfd2d.md, Eqs. 1-2)."""

import scipy.sparse as sp

from .pml import stretch_factors


def second_difference(n, h, k, n_pml):
    """(1/s) d/dx (1/s d/dx) on n nodes with u = 0 just outside: the l_i, r_i of docs/fdfd1d.md, Eq. 2."""
    s, s_half = stretch_factors(n, n_pml, h, k)
    left = 1 / (s * s_half[:-1] * h**2)
    right = 1 / (s * s_half[1:] * h**2)
    return sp.diags([left[1:], -left - right, right[:-1]], [-1, 0, 1], format="csc")


def laplacian_2d(nx, ny, h, k, n_pml):
    """Stretched Laplacian on an (nx, ny) grid, flattened in C order (docs/fdfd2d.md, Eq. 1)."""
    Dxx = second_difference(nx, h, k, n_pml)
    Dyy = second_difference(ny, h, k, n_pml)
    return (sp.kron(Dxx, sp.identity(ny)) + sp.kron(sp.identity(nx), Dyy)).tocsc()