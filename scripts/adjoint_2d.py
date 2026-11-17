"""Figure: forward field, adjoint field and gradient for a mode transmission, plus the 2D V plot."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse.linalg import splu, spsolve

from adjoint_helmholtz.adjoint import intensity_and_gradient
from adjoint_helmholtz.helmholtz2d import helmholtz_matrix
from adjoint_helmholtz.modes import discrete_beta, solve_modes
from adjoint_helmholtz.ports import mode_monitor, mode_source

k0, h, n_pml = 2 * np.pi / 1.55, 0.04, 15
nx, ny = round(3.0 / h) + 1 + 2 * n_pml, round(2.4 / h) + 1 + 2 * n_pml
y = np.arange(ny) * h
eps_r = np.tile(np.where(np.abs(y - y[-1] / 2) <= 0.2, 12.1, 2.1), (nx, 1))
design = (slice(n_pml + 20, nx - n_pml - 20), slice(ny // 2 - 15, ny // 2 + 15))
eps_r[design] = 2.1 + 10 * np.random.default_rng(1).random(eps_r[design].shape)

beta2, e = solve_modes(eps_r[0], h, k0, n_pml)
A = helmholtz_matrix(eps_r, h, k0, n_pml)
b = mode_source(A, eps_r.shape, h, e[:, 0], discrete_beta(beta2[0], h), axis=0, index=n_pml + 5)
c = mode_monitor(eps_r.shape, e[:, 0], axis=0, index=nx - n_pml - 5)

J, grad = intensity_and_gradient(A, b, c, k0)
lu = splu(A)
u = lu.solve(b)
lam = lu.solve(np.conj(np.vdot(c, u)) * np.conj(c), trans="T")


def transmission(eps):
    return abs(np.vdot(c, spsolve(helmholtz_matrix(eps, h, k0, n_pml), b))) ** 2


pixel = (nx // 2, ny // 2)
steps = np.logspace(-10, 0, 31)
errors = []
for step in steps:
    plus, minus = eps_r.copy(), eps_r.copy()
    plus[pixel] += step
    minus[pixel] -= step
    fd = (transmission(plus) - transmission(minus)) / (2 * step)
    errors.append(abs(fd - grad.reshape(nx, ny)[pixel]) / np.max(np.abs(grad)))

fig, ax = plt.subplots(1, 4, figsize=(17, 3.6))
extent = [0, (ny - 1) * h, (nx - 1) * h, 0]
panels = [
    (np.abs(u.reshape(nx, ny)), "$|u|$ forward", "magma"),
    (np.abs(lam.reshape(nx, ny)), "$|\\lambda|$ adjoint", "magma"),
    (grad.reshape(nx, ny), f"$dJ/d\\varepsilon_r$  (J = {J:.3f})", "RdBu"),
]
for axis, (image, title, cmap) in zip(ax, panels):
    limit = np.max(np.abs(image))
    vmin = -limit if cmap == "RdBu" else 0
    im = axis.imshow(image, extent=extent, cmap=cmap, vmin=vmin, vmax=limit)
    axis.set_title(title, fontsize=9)
    axis.set_xlabel("y (um)")
    axis.set_ylabel("x (um)")
    fig.colorbar(im, ax=axis)
ax[3].loglog(steps, errors, "s-")
ax[3].set_xlabel("step in $\\varepsilon_r$")
ax[3].set_ylabel("|FD - adjoint| / max|grad|")
ax[3].set_title(f"central differences at pixel {pixel}", fontsize=9)
fig.tight_layout()
fig.savefig("docs/figures/adjoint_2d.png", dpi=130)
plt.show()