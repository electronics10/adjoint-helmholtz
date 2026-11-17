"""The bend's full chain rho -> filter -> projection -> eps_r -> T, against central differences."""

import numpy as np

from adjoint_helmholtz.bend import Bend


def test_full_chain_gradient_matches_central_differences():
    bend = Bend(h_um=0.05)  # coarse grid: the check is about the chain rule, not accuracy
    rng = np.random.default_rng(0)
    rho = rng.uniform(0.3, 0.7, bend.nd**2)
    step = 1e-4
    for beta in [1, 16]:
        _, grad = bend.objective(rho, beta)
        for p in rng.choice(bend.nd**2, 8, replace=False):
            plus, minus = rho.copy(), rho.copy()
            plus[p] += step
            minus[p] -= step
            fd = (bend.objective(plus, beta)[0] - bend.objective(minus, beta)[0]) / (2 * step)
            assert abs(fd - grad[p]) < 1e-6 * np.max(np.abs(grad)), (beta, p)