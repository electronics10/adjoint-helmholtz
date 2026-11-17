"""Adjoint gradient against JAX automatic differentiation of a dense solve, on a small 2D grid."""

import jax
import jax.numpy as jnp
import numpy as np

from adjoint_helmholtz.adjoint import intensity_and_gradient
from adjoint_helmholtz.helmholtz2d import helmholtz_matrix, point_source

jax.config.update("jax_enable_x64", True)

K0, H, N_PML = 2 * np.pi, 1 / 10, 4
NX, NY = 20, 16


def test_gradient_matches_jax():
    rng = np.random.default_rng(3)
    eps_r = 1 + 2 * rng.random((NX, NY))
    b = point_source((NX, NY), 6, 8, H)
    c = np.zeros(NX * NY, complex)
    c[14 * NY + 5 : 14 * NY + 11] = rng.random(6) + 1j * rng.random(6)

    _, grad = intensity_and_gradient(helmholtz_matrix(eps_r, H, K0, N_PML), b, c, K0)

    # A = L + k0^2 diag(eps_r): take L from our code with eps_r = 0, and let JAX differentiate the rest.
    L = jnp.asarray(helmholtz_matrix(np.zeros((NX, NY)), H, K0, N_PML).toarray())

    def J(eps_flat):
        u = jnp.linalg.solve(L + K0**2 * jnp.diag(eps_flat), b)
        return jnp.abs(jnp.vdot(c, u)) ** 2

    grad_jax = np.asarray(jax.grad(J)(jnp.asarray(eps_r.ravel())))
    assert np.max(np.abs(grad_jax - grad)) < 1e-10 * np.max(np.abs(grad))