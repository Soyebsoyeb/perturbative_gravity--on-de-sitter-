"""Spectra built from numerically evolved modes."""

from __future__ import annotations

import numpy as np

from ..dynamics.numerical_solvers import SolverConfig, solve_mode
from ..dynamics import analytic_solutions as an
from .tensor_power import tensor_power_from_mode, power_exact


def spectrum_at_fixed_eta(ks, eta_end=-1e-2, x_start=100.0, config=None, formulation="conformal",
                          ic_name="bd_exact", bogoliubov=None):
    """Delta_t^2(k) at a common conformal time eta_end for each k in ks.

    Each mode is started at x_start = -k eta_0 (so eta_0 differs per k) and evolved to
    the same final conformal time.  Returns dict(k, numeric, exact, rel_err, x_end).
    """
    cfg = config or SolverConfig(rtol=1e-10, atol=1e-13, n_out=20)
    num, ex = [], []
    for k in ks:
        x_end = -k * eta_end
        if x_end >= x_start:
            raise ValueError(f"mode k={k} has x_end >= x_start; choose smaller k or larger x_start")
        s = solve_mode(k, x_start, x_end, formulation, cfg, ic_name=ic_name, bogoliubov=bogoliubov)
        num.append(tensor_power_from_mode(k, s.gamma[-1]))
        ex.append(power_exact(k, eta_end, cfg.H, cfg.M_pl))
    num, ex = np.array(num), np.array(ex)
    return dict(k=np.asarray(ks, float), numeric=num, exact=ex, rel_err=np.abs(num - ex) / ex,
                x_end=-np.asarray(ks, float) * eta_end)


def amplitude_evolution(k=1.0, x_start=100.0, x_end=1e-3, config=None):
    """|gamma|/gamma_ref versus x for the numerical and exact modes."""
    cfg = config or SolverConfig(n_out=800)
    s = solve_mode(k, x_start, x_end, "conformal", cfg)
    ref = an.gamma_ref(k, cfg.H, cfg.M_pl)
    return dict(x=s.x, numeric=np.abs(s.gamma) / ref,
                exact=np.sqrt(1 + s.x**2), sub=s.x, solution=s)


def horizon_crossing_amplitude(k=1.0, config=None):
    """|gamma| / gamma_ref at x = 1; exact value sqrt(2)."""
    cfg = config or SolverConfig(n_out=5)
    s = solve_mode(k, 50.0, 1.0, "conformal", cfg)
    return float(np.abs(s.gamma[-1]) / an.gamma_ref(k, cfg.H, cfg.M_pl))
