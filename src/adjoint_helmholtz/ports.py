"""One-way mode source and mode-overlap monitor on a grid line (docs/modes.md, Eqs. 4-6).

The line is `index` along `axis` (0: a line of constant x, the mode travels along x; 1: along y).
`direction` is +1 for a mode travelling towards larger index, -1 towards smaller.
"""

import numpy as np


def incident_field(shape, h, profile, beta_d, axis, index, direction):
    """u_inc = e(t) exp(i direction beta_d (s - s_index)) on the whole grid, shape `shape` (Eq. 4)."""
    s = (np.arange(shape[axis]) - index) * h
    phase = np.exp(1j * direction * beta_d * s)
    return np.outer(phase, profile) if axis == 0 else np.outer(profile, phase)


def total_field_mask(shape, axis, index, direction):
    """Q: True on the side of the line the mode travels into, the line itself included."""
    s = np.arange(shape[axis]) - index
    q = direction * s >= 0
    return np.broadcast_to(q[:, None] if axis == 0 else q[None, :], shape)


def mode_source(A, shape, h, profile, beta_d, axis, index, direction=+1):
    """b = (A Q - Q A) u_inc: the field is u_inc inside Q and zero outside (Eq. 5)."""
    u_inc = incident_field(shape, h, profile, beta_d, axis, index, direction).ravel()
    Q = total_field_mask(shape, axis, index, direction).ravel()
    return A @ (Q * u_inc) - Q * (A @ u_inc)


def mode_monitor(shape, profile, axis, index):
    """c such that c^H u is the amplitude of `profile` on the line (Eq. 6)."""
    c = np.zeros(shape, complex)
    line = profile / np.vdot(profile, profile)
    if axis == 0:
        c[index, :] = line
    else:
        c[:, index] = line
    return c.ravel()