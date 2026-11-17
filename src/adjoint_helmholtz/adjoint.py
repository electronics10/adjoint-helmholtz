"""Adjoint gradient of J = |c^H u|^2 with respect to eps_r (docs/adjoint.md)."""

import numpy as np
from scipy.sparse.linalg import splu


def intensity_and_gradient(A, b, c, k0):
    """J = |c^H u|^2 with A u = b, and dJ/d eps_r at every node.

    Assumes dA/d eps_r,i = k0^2 e_i e_i^T (docs/adjoint.md, Eq. 1).
    """
    lu = splu(A.tocsc())  # one factorization for both solves
    u = lu.solve(b)
    a = np.vdot(c, u)  # a = c^H u
    g = np.conj(a) * np.conj(c)  # g = dJ/du, Eq. 3
    lam = lu.solve(g, trans="T")  # A^T lambda = g, Eq. 4
    return abs(a) ** 2, -2 * k0**2 * np.real(lam * u)  # Eq. 5