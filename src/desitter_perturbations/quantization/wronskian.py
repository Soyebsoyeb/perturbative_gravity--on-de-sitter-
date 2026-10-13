"""Wronskians of mode functions."""

from __future__ import annotations

import numpy as np

from ..dynamics import analytic_solutions as an


def canonical_wronskian(v, dv):
    """v v*' - v* v'  (equals i for correctly normalised modes)."""
    return v * np.conj(dv) - np.conj(v) * dv


def gamma_wronskian(gamma, dgamma, z):
    """z^2 (gamma gamma*' - gamma* gamma')  (equals i)."""
    return z**2 * (gamma * np.conj(dgamma) - np.conj(gamma) * dgamma)


def bd_wronskian_deviation(k, eta, H=1.0, M_pl=1.0):
    g = an.gamma_bd(k, eta, H, M_pl)
    dg = an.dgamma_bd(k, eta, H, M_pl)
    W = gamma_wronskian(g, dg, an.z_T(eta, H, M_pl))
    return np.abs(W - 1j)


def conservation_from_solution(sol):
    """Max |W - i| along a numerical solution (should be at the solver-error level)."""
    z = an.z_T(sol.eta, sol.config.H, sol.config.M_pl)
    W = gamma_wronskian(sol.gamma, sol.dgamma_deta, z)
    return float(np.max(np.abs(W - 1j)))
