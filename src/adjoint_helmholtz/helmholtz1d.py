"""1D Helmholtz operator with PML: A u = b discretizes u'' + k0^2 eps_r u = b (docs/fdfd1d.md)."""

import numpy as np
import scipy.sparse as sp

from .pml import stretch_factors


def helmholtz_matrix(eps_r, h, k0, n_pml):
    """Sparse A, shape (n, n), with the PML rows of docs/fdfd1d.md, Eq. 2. Rows are not scaled."""
    n = len(eps_r)
    s, s_half = stretch_factors(n, n_pml, h, k0)
    left = 1 / (s * s_half[:-1] * h**2)  # coefficient of u_{i-1}
    right = 1 / (s * s_half[1:] * h**2)  # coefficient of u_{i+1}
    diag = k0**2 * eps_r - left - right
    return sp.diags([left[1:], diag, right[:-1]], [-1, 0, 1], format="csc")


def point_source(n, j, h):
    """Unit point source -delta at node j: b_j = -1/h (conventions.md, item 2)."""
    b = np.zeros(n, complex)
    b[j] = -1 / h
    return b