"""Density filter and projection for topology optimization (docs/topopt.md).

rho (design variables in [0, 1]) -> rho_t = W rho (filtered) -> rho_p = P(rho_t) (projected)
"""

import numpy as np
import scipy.sparse as sp


def conic_filter(shape, radius_px):
    """Sparse W, row-normalized, with weights max(0, R - distance) over the design grid (Eq. 1).

    shape is the design grid (nx, ny), flattened in C order. Near the edges a row averages over
    fewer pixels, because there are no pixels outside the design region.
    """
    nx, ny = shape
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    reach = int(np.ceil(radius_px))
    rows, cols, weights = [], [], []
    for di in range(-reach, reach + 1):
        for dj in range(-reach, reach + 1):
            weight = radius_px - np.hypot(di, dj)
            if weight <= 0:
                continue
            inside = (i + di >= 0) & (i + di < nx) & (j + dj >= 0) & (j + dj < ny)
            rows.append((i * ny + j)[inside])
            cols.append(((i + di) * ny + (j + dj))[inside])
            weights.append(np.full(inside.sum(), weight))
    n = nx * ny
    W = sp.csr_matrix((np.concatenate(weights), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))
    return sp.diags(1 / np.asarray(W.sum(axis=1)).ravel()) @ W


def project(rho_t, beta, eta=0.5):
    """Smoothed step at eta, with steepness beta; maps 0 to 0 and 1 to 1 (Eq. 2)."""
    return (np.tanh(beta * eta) + np.tanh(beta * (rho_t - eta))) / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))


def project_derivative(rho_t, beta, eta=0.5):
    """d project / d rho_t (Eq. 3)."""
    return beta * (1 - np.tanh(beta * (rho_t - eta)) ** 2) / (np.tanh(beta * eta) + np.tanh(beta * (1 - eta)))


def binary_fraction(rho_p, tolerance=0.05):
    """Fraction of pixels within tolerance of 0 or 1."""
    return float(np.mean((rho_p <= tolerance) | (rho_p >= 1 - tolerance)))