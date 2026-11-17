"""Optimize the 90-degree bend with beta continuation and L-BFGS-B (docs/bend.md).

Writes results/bend.npz (not in git) and docs/figures/bend.png.
"""

import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from adjoint_helmholtz.bend import Bend
from adjoint_helmholtz.topopt import binary_fraction, project

BETAS = [1, 2, 4, 8, 16, 32, 64, 128]
ITERATIONS_PER_BETA = 40

bend = Bend()
T_circle, _ = bend.transmission(bend.quarter_circle)
print(f"grid {bend.shape}, design {bend.nd} x {bend.nd}, quarter circle T = {T_circle:.4f}")

rho = np.full(bend.nd**2, 0.5)
history = []  # (beta, T) at every objective evaluation
for beta in BETAS:
    start = time.time()

    def negative_T(r):
        T, grad = bend.objective(r, beta)
        history.append((beta, T))
        return -T, -grad

    result = minimize(negative_T, rho, jac=True, method="L-BFGS-B", bounds=[(0, 1)] * rho.size,
                      options={"maxiter": ITERATIONS_PER_BETA})
    rho = result.x
    rho_p = project(bend.W @ rho, beta)
    print(f"beta {beta:3d}: T = {-result.fun:.4f}, binary fraction = {binary_fraction(rho_p):.3f}, "
          f"{result.nfev} evaluations, {time.time() - start:.0f} s")

rho_p = project(bend.W @ rho, BETAS[-1])
T_continuous, _ = bend.transmission(rho_p)
T_binary, u = bend.transmission((rho_p > 0.5).astype(float))
print(f"final: T = {T_continuous:.4f} (continuous), {T_binary:.4f} (thresholded at 0.5), "
      f"binary fraction {binary_fraction(rho_p):.3f}; quarter circle {T_circle:.4f}")
np.savez("results/bend.npz", rho=rho, rho_p=rho_p, history=np.array(history))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
history = np.array(history)
ax[0].plot(history[:, 1])
for beta in BETAS:
    first = np.argmax(history[:, 0] == beta)
    ax[0].axvline(first, color="gray", lw=0.5)
    ax[0].text(first, 0.03, f"$\\beta$={beta}", fontsize=7)
ax[0].axhline(T_circle, color="tab:red", ls="--", label=f"quarter circle, T = {T_circle:.3f}")
ax[0].set_xlabel("objective evaluation")
ax[0].set_ylabel("T")
ax[0].set_ylim(0, 1.02)
ax[0].legend(loc="lower right")
ax[0].set_title("convergence over the beta steps")
eps_r = bend.eps_r((rho_p > 0.5).astype(float))
extent = [0, bend.shape[0] * bend.h, 0, bend.shape[1] * bend.h]
ax[1].imshow(eps_r.T, origin="lower", cmap="Greys", extent=extent)
ax[1].set_title(f"thresholded design, T = {T_binary:.4f}")
ax[2].imshow(np.abs(u.reshape(bend.shape)).T, origin="lower", cmap="magma", extent=extent)
ax[2].contour(eps_r.T, levels=[7], colors="w", linewidths=0.5, extent=extent)
ax[2].set_title("$|u|$, thresholded design")
for a in ax[1:]:
    a.set_xlabel("x (um)")
    a.set_ylabel("y (um)")
fig.tight_layout()
fig.savefig("docs/figures/bend.png", dpi=120)
plt.show()