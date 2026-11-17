"""Modes of a straight waveguide, from the eps_r on one grid line (docs/modes.md, Eqs. 1-3)."""

import numpy as np

from .operators import second_difference


def solve_modes(eps_line, h, k0, n_pml, num_modes=1):
    """Guided modes, largest beta^2 first. Returns beta2, shape (num_modes,), and profiles, shape (n, num_modes).

    Each profile e has sum |e|^2 h = 1 and a real positive peak.
    """
    D = second_difference(len(eps_line), h, k0, n_pml).toarray()
    beta2, e = np.linalg.eig(D + np.diag(k0**2 * eps_line))  # Eq. 2
    guided = np.abs(beta2.imag) < 1e-6 * np.abs(beta2.real)  # drops modes that live in the PML
    beta2, e = beta2[guided], e[:, guided]
    order = np.argsort(-beta2.real)[:num_modes]
    beta2, e = beta2[order], e[:, order]
    for m in range(num_modes):
        peak = e[np.argmax(np.abs(e[:, m])), m]
        e[:, m] *= (np.conj(peak) / abs(peak)) / np.sqrt(np.sum(np.abs(e[:, m]) ** 2) * h)
    return beta2, e


def discrete_beta(beta2, h):
    """Grid wavenumber along the guide, so that the mode solves the discrete 2D equation (Eq. 3)."""
    return np.arccos(1 - beta2 * h**2 / 2) / h