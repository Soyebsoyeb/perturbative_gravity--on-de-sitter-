"""Convergence and stability studies for the tensor-mode solvers."""

from __future__ import annotations

import numpy as np

from . import analytic_solutions as an
from . import initial_conditions as ic
from .numerical_solvers import SolverConfig, solve_mode, solve_rk4_efold
from .residuals import error_vs_exact


def observed_order(h, err):
    """Least-squares slope of log(err) vs log(h)."""
    h = np.asarray(h, float)
    err = np.asarray(err, float)
    m = err > 0
    return float(np.polyfit(np.log(h[m]), np.log(err[m]), 1)[0])


def tolerance_study(k=1.0, rtols=(1e-4, 1e-6, 1e-8, 1e-10, 1e-12), x_start=100.0, x_end=1e-3,
                    formulation="conformal", method="DOP853"):
    out = []
    for rtol in rtols:
        cfg = SolverConfig(method=method, rtol=rtol, atol=rtol * 1e-3, n_out=600)
        s = solve_mode(k, x_start, x_end, formulation, cfg)
        e = error_vs_exact(s)
        out.append(dict(rtol=rtol, max_err=float(np.max(e)), final_err=float(e[-1]), nfev=s.nfev))
    return out


def error_trajectories(k=1.0, rtols=(1e-4, 1e-7, 1e-10), x_start=100.0, x_end=1e-3,
                       formulation="conformal", method="DOP853"):
    """Relative error along the trajectory (shows behaviour at horizon crossing x = 1)."""
    res = {}
    for rtol in rtols:
        cfg = SolverConfig(method=method, rtol=rtol, atol=rtol * 1e-3, n_out=800)
        s = solve_mode(k, x_start, x_end, formulation, cfg)
        res[rtol] = (s.x, error_vs_exact(s))
    return res


def method_comparison(k=1.0, rtol=1e-8, methods=("RK45", "RK23", "DOP853", "Radau", "BDF", "LSODA"),
                      x_start=100.0, x_end=1e-3):
    out = []
    for m in methods:
        cfg = SolverConfig(method=m, rtol=rtol, atol=rtol * 1e-3, n_out=600)
        # LSODA/BDF/Radau accept real systems; our state is real
        s = solve_mode(k, x_start, x_end, "conformal", cfg)
        e = error_vs_exact(s)
        out.append(dict(method=m, max_err=float(np.max(e)), nfev=s.nfev, success=s.success))
    return out


def formulation_comparison(k=1.0, rtol=1e-10, x_start=100.0, x_end=1e-3):
    cfg = SolverConfig(rtol=rtol, atol=rtol * 1e-3, n_out=600)
    sols = {f: solve_mode(k, x_start, x_end, f, cfg) for f in ("conformal", "cosmic", "efold", "canonical")}
    errs = {f: float(np.max(error_vs_exact(s))) for f, s in sols.items()}
    return sols, errs


def rk4_order_study(k=1.0, x_start=30.0, x_end=1e-2, steps=(250, 500, 1000, 2000, 4000)):
    errs, hs = [], []
    ex = an.gamma_bd(k, -x_end / k)
    for n in steps:
        _, g, _ = solve_rk4_efold(k, x_start, x_end, n)
        errs.append(abs(g - ex) / abs(ex))
        hs.append(np.log(x_start / x_end) / n)
    return np.array(hs), np.array(errs), observed_order(hs, errs)


def start_x_study(k=1.0, x_starts=(5, 10, 20, 50, 100, 200, 500, 1000), x_end=1e-3, rtol=1e-11):
    """Effect of imposing only the leading positive-frequency asymptotics at finite x_start.

    Returns rows with the numerical power ratio P/P_BD at x_end, the exact Bogoliubov
    prediction |alpha - beta e^{...}|^2 evaluated from the exact BD modes, and |beta|.
    """
    rows = []
    for xs in x_starts:
        cfg = SolverConfig(rtol=rtol, atol=rtol * 1e-3, n_out=50)
        s = solve_mode(k, xs, x_end, "conformal", cfg, ic_name="bd_asymptotic")
        ex = an.gamma_bd(k, s.eta[-1])
        ratio = abs(s.gamma[-1]) ** 2 / abs(ex) ** 2
        eta0 = -xs / k
        al, be = ic.bogoliubov_coefficients(k, eta0, *ic.bd_asymptotic(k, eta0))
        # exact prediction: gamma(eta) = al u + be u*, ratio = |al u + be u*|^2 / |u|^2
        u = ex
        pred = abs(al * u + be * np.conj(u)) ** 2 / abs(u) ** 2
        rows.append(dict(x_start=xs, ratio=float(ratio), predicted=float(pred), abs_beta=float(abs(be))))
    return rows
