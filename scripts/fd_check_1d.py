"""V plot: relative error of finite differences against the adjoint gradient, versus step size."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import spsolve

from adjoint_helmholtz.adjoint import intensity_and_gradient
from adjoint_helmholtz.helmholtz1d import helmholtz_matrix, point_source

k0, h, n_pml, n = 2 * np.pi, 1 / 20, 10, 101
eps_r = np.ones(n)
eps_r[40:60] = 2.0
b = point_source(n, 20, h)
c = np.zeros(n)
c[80] = 1
i = 45  # node to check


def intensity(eps):
    u = spsolve(helmholtz_matrix(eps, h, k0, n_pml), b)
    return abs(np.vdot(c, u)) ** 2


_, grad = intensity_and_gradient(helmholtz_matrix(eps_r, h, k0, n_pml), b, c, k0)
steps = np.logspace(-12, -1, 45)
forward, central = [], []
for step in steps:
    plus, minus = eps_r.copy(), eps_r.copy()
    plus[i] += step
    minus[i] -= step
    J0, Jp, Jm = intensity(eps_r), intensity(plus), intensity(minus)
    forward.append(abs((Jp - J0) / step - grad[i]) / abs(grad[i]))
    central.append(abs((Jp - Jm) / (2 * step) - grad[i]) / abs(grad[i]))

plt.loglog(steps, forward, "o-", label="forward (error ~ step)")
plt.loglog(steps, central, "s-", label="central (error ~ step$^2$)")
plt.xlabel("step in $\\varepsilon_r$")
plt.ylabel("relative error vs adjoint")
plt.legend()
plt.title(f"Finite differences vs adjoint gradient, node {i}")
plt.savefig("docs/figures/fd_check_1d.png", dpi=150)
plt.show()