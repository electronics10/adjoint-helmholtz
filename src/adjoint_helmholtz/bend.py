"""A 90-degree waveguide bend as a topology-optimization problem (docs/bend.md).

Si (eps 12.1) in SiO2 (eps 2.1) at 1.55 um. The input guide enters a square design region from the
left (along +x); the output guide leaves it from the bottom (along -y). T is the power in the output
guide's fundamental mode, for a launched amplitude of 1.
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from .adjoint import intensity_and_gradient
from .helmholtz2d import helmholtz_matrix
from .modes import discrete_beta, solve_modes
from .ports import mode_monitor, mode_source
from .topopt import conic_filter, project, project_derivative

WAVELENGTH_UM = 1.55
EPS_MIN, EPS_MAX = 2.1, 12.1
K0 = 2 * np.pi / WAVELENGTH_UM


class Bend:
    def __init__(self, h_um=0.02, width_um=0.4, design_um=2.4, pml_um=0.8, lead_um=1.0, gap_um=0.5, filter_um=0.1):
        self.h = h_um
        self.n_pml = round(pml_um / h_um)
        n = round((2 * pml_um + lead_um + design_um + gap_um) / h_um) + 1
        self.shape = (n, n)
        start = round((pml_um + lead_um) / h_um)  # first design index, the same in x and y
        self.nd = round(design_um / h_um)
        self.design = (slice(start, start + self.nd), slice(start, start + self.nd))

        # Background: two straight guides meeting the design region. Guide axes cross at its centre.
        x = np.arange(n) * h_um
        X, Y = np.meshgrid(x, x, indexing="ij")
        axis_um = (start + self.nd / 2) * h_um
        edge_um = start * h_um
        self.eps_background = np.full(self.shape, EPS_MIN)
        self.eps_background[(np.abs(Y - axis_um) <= width_um / 2) & (X < edge_um)] = EPS_MAX
        self.eps_background[(np.abs(X - axis_um) <= width_um / 2) & (Y < edge_um)] = EPS_MAX

        # A without the k0^2 eps_r term, built once (docs/fdfd2d.md, Eq. 2).
        self.L = helmholtz_matrix(np.zeros(self.shape), h_um, K0, self.n_pml)

        i_source = self.n_pml + 5
        beta2, e = solve_modes(self.eps_background[i_source, :], h_um, K0, self.n_pml)
        A = self.matrix(self.eps_background)
        self.b = mode_source(A, self.shape, h_um, e[:, 0], discrete_beta(beta2[0], h_um), axis=0, index=i_source)
        j_monitor = self.n_pml + 10
        beta2, e = solve_modes(self.eps_background[:, j_monitor], h_um, K0, self.n_pml)
        self.c = mode_monitor(self.shape, e[:, 0], axis=1, index=j_monitor)

        self.W = conic_filter((self.nd, self.nd), filter_um / h_um)

        # Baseline: a quarter-circle guide centred on the design region's bottom-left corner.
        radius_um = axis_um - edge_um
        ring = np.abs(np.hypot(X - edge_um, Y - edge_um) - radius_um) <= width_um / 2
        self.quarter_circle = ring[self.design].astype(float)

    def matrix(self, eps_r):
        return (self.L + sp.diags(K0**2 * eps_r.ravel())).tocsc()

    def eps_r(self, rho_p):
        """Background with the design region filled by eps_min + rho_p (eps_max - eps_min) (docs/topopt.md, Eq. 4)."""
        eps_r = self.eps_background.copy()
        eps_r[self.design] = EPS_MIN + rho_p.reshape(self.nd, self.nd) * (EPS_MAX - EPS_MIN)
        return eps_r

    def transmission(self, rho_p):
        """T and the field u for a projected density rho_p on the design grid."""
        u = spsolve(self.matrix(self.eps_r(rho_p)), self.b)
        return abs(np.vdot(self.c, u)) ** 2, u

    def objective(self, rho, beta):
        """T and dT/drho for the full chain rho -> filter -> projection -> eps_r -> T (docs/topopt.md, Eq. 5)."""
        rho_t = self.W @ rho
        rho_p = project(rho_t, beta)
        T, dT_deps = intensity_and_gradient(self.matrix(self.eps_r(rho_p)), self.b, self.c, K0)
        dT_deps_design = dT_deps.reshape(self.shape)[self.design].ravel()
        return T, self.W.T @ (project_derivative(rho_t, beta) * (EPS_MAX - EPS_MIN) * dT_deps_design)