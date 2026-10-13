"""Diagnostics that do not rely on knowing the exact solution.

* ``equation_residuals``: splits the second-order equation into two first-order
  consistency checks evaluated with an independent cubic-spline derivative of the
  stored output (kinematic: d gamma/d eta = gamma'; dynamic: d gamma'/d eta = RHS).
* ``wronskian_deviation``: relative drift of the conserved Wronskian
  z^2 (gamma gamma*' - gamma* gamma') = i  from its exact value.
* ``cross_formulation``: pointwise relative differences between formulations.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from . import analytic_solutions as an


def equation_residuals(sol):
    """Return (kinematic, dynamic) residual arrays, normalised by local scales."""
    eta, g, dg, k = sol.eta, sol.gamma, sol.dgamma_deta, sol.k
    # eta is increasing along the solution (negative -> 0^-)
    cs_g = CubicSpline(eta, g)
    cs_dg = CubicSpline(eta, dg)
    kin = cs_g(eta, 1) - dg
    dyn = cs_dg(eta, 1) - ((2.0 / eta) * dg - k * k * g)
    scale = np.maximum(np.abs(g) * k, np.abs(dg)) + 1e-300
    scale2 = np.maximum(np.abs(g) * k * k, np.abs(dg) * 2.0 / np.abs(eta)) + 1e-300
    return np.abs(kin) / scale, np.abs(dyn) / scale2


def wronskian_deviation(sol):
    """Relative deviation of gamma gamma*' - gamma* gamma' from 2 i H^2 eta^2 / M^2."""
    H, M = sol.config.H, sol.config.M_pl
    W = sol.gamma * np.conj(sol.dgamma_deta) - np.conj(sol.gamma) * sol.dgamma_deta
    W_exact = an.wronskian_gamma_exact(sol.eta, H, M)
    return np.abs(W - W_exact) / np.abs(W_exact)


def error_vs_exact(sol):
    """Pointwise relative error |gamma - gamma_BD| / |gamma_BD| (exact BD initial data only)."""
    ex = an.gamma_bd(sol.k, sol.eta, sol.config.H, sol.config.M_pl)
    return np.abs(sol.gamma - ex) / np.abs(ex)


def cross_formulation(solutions, reference="conformal"):
    """Dict formulation -> max over x of |gamma_f - gamma_ref| / |gamma_ref|."""
    ref = solutions[reference].gamma
    return {name: float(np.max(np.abs(s.gamma - ref) / np.abs(ref)))
            for name, s in solutions.items() if name != reference}
