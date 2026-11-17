"""2D Helmholtz operator with PML: A u = b discretizes (nabla^2 + k0^2 eps_r) u = b (docs/fdfd2d.md)."""

import numpy as np
import scipy.sparse as sp

from .operators import laplacian_2d


def helmholtz_matrix(eps_r, h, k0, n_pml):
    """Sparse A for eps_r of shape (nx, ny); u is flattened in C order (docs/fdfd2d.md, Eq. 2)."""
    nx, ny = eps_r.shape
    return (laplacian_2d(nx, ny, h, k0, n_pml) + sp.diags(k0**2 * eps_r.ravel())).tocsc()


def point_source(shape, i, j, h):
    """Unit point source -delta at node (i, j): b = -1/h^2 there (conventions.md, item 2)."""
    b = np.zeros(shape, complex)
    b[i, j] = -1 / h**2
    return b.ravel()