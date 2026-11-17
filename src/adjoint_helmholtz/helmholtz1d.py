"""1D Helmholtz operator with PML: A u = b discretizes u'' + k0^2 eps_r u = b (docs/fdfd1d.md)."""

import numpy as np
import scipy.sparse as sp

from .operators import second_difference


def helmholtz_matrix(eps_r, h, k0, n_pml):
    """Sparse A, shape (n, n), with the PML rows of docs/fdfd1d.md, Eq. 2. Rows are not scaled."""
    return (second_difference(len(eps_r), h, k0, n_pml) + sp.diags(k0**2 * eps_r)).tocsc()


def point_source(n, j, h):
    """Unit point source -delta at node j: b_j = -1/h (conventions.md, item 2)."""
    b = np.zeros(n, complex)
    b[j] = -1 / h
    return b