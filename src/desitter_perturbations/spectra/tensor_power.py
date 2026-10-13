"""Dimensionless tensor power spectrum.

Definition (declared convention):

    <gamma_ij(x) gamma_ij(x)> = int dln k  Delta_t^2(k),
    gamma_ij(x) = sum_lam int d^3k/(2 pi)^3 e^lam_ij(khat) gamma_lam(k) e^{i k.x}
                  with <gamma_lam(k) gamma_lam'(k')> = (2 pi)^3 delta^3(k + k') delta_{lam lam'} |gamma_k|^2.

With n_pol polarizations of normalisation e^lam_ij e^lam_ij = N_e (=2):

    Delta_t^2(k) = n_pol N_e k^3 |gamma_k|^2 / (2 pi^2)  =  2 k^3 |gamma_k|^2 / pi^2   (n_pol = N_e = 2).

For the exact BD mode, |gamma_k|^2 = H^2 (1 + k^2 eta^2) / (M_pl^2 k^3), so

    Delta_t^2(k, eta) = (2 H^2 / (pi^2 M_pl^2)) (1 + k^2 eta^2)  ->  2 H^2/(pi^2 M_pl^2)  (superhorizon).

The check is on the *implementation*: this formula is verified numerically against
the numerically evolved modes (``mode_spectrum``) and the prefactor follows from the
same polarization convention that fixes z_T (see docs/power_spectrum.md).
"""

from __future__ import annotations

import numpy as np

from ..conventions import POLARIZATION_NORM


def tensor_power_from_mode(k, gamma_k, n_pol: int = 2, pol_norm: float = POLARIZATION_NORM):
    k = np.asarray(k, dtype=float)
    return n_pol * pol_norm * k**3 * np.abs(gamma_k) ** 2 / (2.0 * np.pi**2)


def power_exact(k, eta, H=1.0, M_pl=1.0):
    k = np.asarray(k, dtype=float)
    eta = np.asarray(eta, dtype=float)
    return 2.0 * H**2 / (np.pi**2 * M_pl**2) * (1.0 + (k * eta) ** 2)


def power_superhorizon(H=1.0, M_pl=1.0):
    return 2.0 * H**2 / (np.pi**2 * M_pl**2)


def power_bogoliubov_superhorizon(alpha, beta, H=1.0, M_pl=1.0):
    """Leading superhorizon spectrum for gamma = alpha u + beta u*: |alpha - beta|^2 times BD."""
    return abs(alpha - beta) ** 2 * power_superhorizon(H, M_pl)


def tensor_tilt(k, P):
    """n_t = d ln P / d ln k by centred differences on a log grid."""
    lk, lp = np.log(k), np.log(P)
    return np.gradient(lp, lk)


def variance_from_spectrum(k, P):
    """int dln k P(k) by the trapezoid rule on the supplied (log-spaced) grid."""
    lk = np.log(k)
    return float(np.trapezoid(P, lk))
