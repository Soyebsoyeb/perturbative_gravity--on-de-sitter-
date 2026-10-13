"""Initial data for tensor modes at x_start = -k eta_0.

* ``bd_exact``       exact Bunch-Davies solution evaluated at eta_0.
* ``bd_asymptotic``  positive-frequency WKB data  v = e^{-i k eta0}/sqrt(2k),
                     v' = -i k v   (what one would impose knowing only the
                     short-wavelength limit), converted to (gamma, gamma').
* ``bogoliubov``     exact modes mixed with Bogoliubov coefficients (alpha, beta),
                     |alpha|^2 - |beta|^2 = 1.

All return (gamma, dgamma/deta) in conformal time.  ``bogoliubov_coefficients``
projects arbitrary data onto the exact BD modes via the Wronskian.
"""

from __future__ import annotations

import numpy as np

from . import analytic_solutions as an


def bd_exact(k, eta0, H=1.0, M_pl=1.0):
    return an.gamma_bd(k, eta0, H, M_pl), an.dgamma_bd(k, eta0, H, M_pl)


def bd_asymptotic(k, eta0, H=1.0, M_pl=1.0):
    k = float(k)
    v = np.exp(-1j * k * eta0) / np.sqrt(2 * k)
    dv = -1j * k * v
    z = an.z_T(eta0, H, M_pl)
    calH = -1.0 / eta0
    g = v / z
    dg = (dv - calH * v) / z
    return g, dg


def bogoliubov(k, eta0, alpha, beta, H=1.0, M_pl=1.0):
    g, dg = bd_exact(k, eta0, H, M_pl)
    return alpha * g + beta * np.conj(g), alpha * dg + beta * np.conj(dg)


def bogoliubov_coefficients(k, eta0, g0, dg0, H=1.0, M_pl=1.0):
    """Project data (g0, dg0) at eta0 onto u = gamma_bd and u*:  g = alpha u + beta u*.

    Uses the Wronskian W[u, u*] = u u*' - u* u' = 2i/(M^2 a^2) =: w:
        alpha =  (g u*' - g' u*) / w,       beta = -(g u' - g' u) / w.
    """
    u, du = bd_exact(k, eta0, H, M_pl)
    w = an.wronskian_gamma_exact(eta0, H, M_pl)
    alpha = (g0 * np.conj(du) - dg0 * np.conj(u)) / w
    beta = -(g0 * du - dg0 * u) / w
    return alpha, beta
