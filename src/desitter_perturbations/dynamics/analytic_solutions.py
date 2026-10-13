"""Exact and asymptotic tensor-mode solutions in flat de Sitter space.

Equation (conformal time, a = -1/(H eta)):   gamma'' - (2/eta) gamma' + k^2 gamma = 0.

Exact Bunch-Davies solution (derived in docs/tensor_dynamics.md):

    gamma_k(eta) = (i H / (M_pl k^{3/2})) (1 + i k eta) exp(-i k eta),
    gamma_k'(eta) = (i H / (M_pl k^{3/2})) k^2 eta exp(-i k eta).

Canonical variable v = z gamma with z = a M_pl / sqrt(2):

    v_k = exp(-i k eta) / sqrt(2k) * (1 - i/(k eta)).

Asymptotics
-----------
subhorizon (x = -k eta >> 1):    gamma -> exp(-i k eta) / (a M_pl sqrt(k))      (|gamma| decays as 1/a)
superhorizon (x << 1):           gamma -> i H / (M_pl k^{3/2})                   (frozen)
"""

from __future__ import annotations

import numpy as np


def gamma_ref(k, H=1.0, M_pl=1.0):
    """Reference amplitude H / (M_pl k^{3/2}) (the exact superhorizon |gamma|)."""
    return H / (M_pl * np.asarray(k, dtype=float) ** 1.5)


def gamma_bd(k, eta, H=1.0, M_pl=1.0):
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    return 1j * gamma_ref(k, H, M_pl) * (1 + 1j * k * eta) * np.exp(-1j * k * eta)


def dgamma_bd(k, eta, H=1.0, M_pl=1.0):
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    return 1j * gamma_ref(k, H, M_pl) * k**2 * eta * np.exp(-1j * k * eta)


def v_bd(k, eta):
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    return np.exp(-1j * k * eta) / np.sqrt(2 * k) * (1 - 1j / (k * eta))


def dv_bd(k, eta):
    """d v/d eta of the exact BD canonical mode."""
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    pref = np.exp(-1j * k * eta) / np.sqrt(2 * k)
    return pref * (-1j * k * (1 - 1j / (k * eta)) + 1j / (k * eta**2))


def z_T(eta, H=1.0, M_pl=1.0):
    """Canonical normalisation factor z_T = a M_pl / sqrt(2) (derived in quantization.canonical_variables)."""
    return (-1.0 / (H * np.asarray(eta, dtype=float))) * M_pl / np.sqrt(2.0)


def gamma_subhorizon(k, eta, H=1.0, M_pl=1.0):
    a = -1.0 / (H * np.asarray(eta, dtype=float))
    return np.exp(-1j * k * np.asarray(eta, dtype=float)) / (a * M_pl * np.sqrt(k))


def gamma_superhorizon(k, H=1.0, M_pl=1.0):
    return 1j * gamma_ref(k, H, M_pl) * np.ones_like(np.asarray(k, dtype=float), dtype=complex)


def abs_gamma_squared_exact(k, eta, H=1.0, M_pl=1.0):
    """|gamma_k|^2 = H^2 (1 + k^2 eta^2) / (M_pl^2 k^3)."""
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    return H**2 * (1 + (k * eta) ** 2) / (M_pl**2 * k**3)


def wronskian_gamma_exact(eta, H=1.0, M_pl=1.0):
    """gamma gamma*' - gamma* gamma' = 2 i / (M_pl^2 a^2) = 2 i H^2 eta^2 / M_pl^2."""
    return 2j * H**2 * np.asarray(eta, dtype=float) ** 2 / M_pl**2


def residual_exact_conformal(k, eta, H=1.0, M_pl=1.0):
    """Residual gamma'' + 2 calH gamma' + k^2 gamma of the exact solution, using the
    analytic second derivative gamma'' = (iH/(M k^{3/2})) k^2 (1 - i k eta) e^{-i k eta}."""
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    g = gamma_bd(k, eta, H, M_pl)
    dg = dgamma_bd(k, eta, H, M_pl)
    d2g = 1j * gamma_ref(k, H, M_pl) * k**2 * (1 - 1j * k * eta) * np.exp(-1j * k * eta)
    return d2g - 2.0 / eta * dg + k**2 * g
