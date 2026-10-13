"""Right-hand sides of the tensor-mode equation in four formulations.

State vectors are real: [Re g, Im g, Re g_dot, Im g_dot] where the dot is the
derivative with respect to the independent variable of that formulation.

conformal   eta :  g'' = (2/eta) g' - k^2 g
cosmic      t   :  g_tt = -3 H g_t - (k^2/a^2) g,       a = exp(H t)
e-fold      N   :  g_NN = -3 g_N - (k^2 e^{-2N}/H^2) g,  N = ln a
canonical   eta :  v'' = -(k^2 - 2/eta^2) v,            gamma = v / z_T

All use H and the wavenumber k in the same (comoving, a(t=0)=1) units.
"""

from __future__ import annotations

import numpy as np

FORMULATIONS = ("conformal", "cosmic", "efold", "canonical")


def rhs_conformal(eta, y, k):
    g = y[0] + 1j * y[1]
    dg = y[2] + 1j * y[3]
    d2g = (2.0 / eta) * dg - k * k * g
    return np.array([dg.real, dg.imag, d2g.real, d2g.imag])


def rhs_cosmic(t, y, k, H):
    g = y[0] + 1j * y[1]
    dg = y[2] + 1j * y[3]
    a2 = np.exp(2.0 * H * t)
    d2g = -3.0 * H * dg - (k * k / a2) * g
    return np.array([dg.real, dg.imag, d2g.real, d2g.imag])


def rhs_efold(N, y, k, H):
    g = y[0] + 1j * y[1]
    dg = y[2] + 1j * y[3]
    d2g = -3.0 * dg - (k * k * np.exp(-2.0 * N) / H**2) * g
    return np.array([dg.real, dg.imag, d2g.real, d2g.imag])


def rhs_canonical(eta, y, k):
    v = y[0] + 1j * y[1]
    dv = y[2] + 1j * y[3]
    d2v = -(k * k - 2.0 / eta**2) * v
    return np.array([dv.real, dv.imag, d2v.real, d2v.imag])


def pack(g, dg):
    return np.array([np.real(g), np.imag(g), np.real(dg), np.imag(dg)], dtype=float)


def unpack(y):
    y = np.asarray(y)
    return y[0] + 1j * y[1], y[2] + 1j * y[3]
