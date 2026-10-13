"""Vacuum states of the free tensor field on de Sitter space.

* Bunch-Davies (BD): v_k -> e^{-i k eta}/sqrt(2k) as k eta -> -infinity.  It is the
  state annihilated by a_k defined with the positive-frequency modes of the
  early-time (subhorizon) Minkowski-like regime.
* Bogoliubov-rotated states: mode = alpha u + beta u*, |alpha|^2 - |beta|^2 = 1.
  The number of BD quanta is |beta|^2; late-time power is multiplied by |alpha - beta|^2.
* Finite-time WKB state: positive-frequency data imposed at finite x_start; equals a
  Bogoliubov state with |beta| = 1/(2 x_start^2) (derived and verified numerically).

What this does *not* establish: that the BD state is the physically correct state of
the real universe, nor anything about interacting (non-free) dynamics.
"""

from __future__ import annotations

import numpy as np

from ..dynamics import analytic_solutions as an
from ..dynamics import initial_conditions as ic


def bd_mode(k, eta, H=1.0, M_pl=1.0):
    return an.gamma_bd(k, eta, H, M_pl)


def bogoliubov_mode(k, eta, alpha, beta, H=1.0, M_pl=1.0):
    u = an.gamma_bd(k, eta, H, M_pl)
    return alpha * u + beta * np.conj(u)


def check_bogoliubov(alpha, beta, tol=1e-12):
    return abs(abs(alpha) ** 2 - abs(beta) ** 2 - 1.0) < tol


def particle_number(beta):
    return abs(beta) ** 2


def wkb_state_beta(x_start):
    """Bogoliubov coefficient magnitude of the leading-WKB initial state at finite x_start."""
    return 1.0 / (2.0 * x_start**2)


def superhorizon_power_factor(alpha, beta):
    """|alpha - beta|^2 multiplying the BD superhorizon spectrum."""
    return abs(alpha - beta) ** 2


def mode_two_point_k(k, eta, alpha=1.0, beta=0.0, H=1.0, M_pl=1.0):
    """|gamma_k(eta)|^2, the (2 pi)^3 delta-stripped equal-time two-point function."""
    return abs(bogoliubov_mode(k, eta, alpha, beta, H, M_pl)) ** 2
