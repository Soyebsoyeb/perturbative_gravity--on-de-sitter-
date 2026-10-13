"""Parameter sweeps: results are returned as dicts of arrays and can be written to CSV."""

from __future__ import annotations

import csv
import os
import numpy as np

from ..dynamics.numerical_solvers import SolverConfig, solve_mode
from ..dynamics import analytic_solutions as an, convergence as cv, initial_conditions as ic
from .mode_spectrum import spectrum_at_fixed_eta
from .tensor_power import tensor_power_from_mode, power_exact, power_superhorizon


def write_csv(path, table: dict):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    keys = list(table.keys())
    n = len(next(iter(table.values())))
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(n):
            w.writerow([f"{np.asarray(table[k])[i]:.12e}" if np.issubdtype(np.asarray(table[k]).dtype, np.number)
                        else table[k][i] for k in keys])


def sweep_wavenumber(n=25, k_min=0.1, k_max=100.0, eta_end=-1e-2, x_start=200.0):
    ks = np.geomspace(k_min, k_max, n)
    return spectrum_at_fixed_eta(ks, eta_end=eta_end, x_start=x_start)


def sweep_late_time(x_ends=np.geomspace(1e-4, 3.0, 20), k=1.0, x_start=100.0):
    cfg = SolverConfig(rtol=1e-11, atol=1e-14, n_out=10)
    num, ex = [], []
    for xe in x_ends:
        s = solve_mode(k, x_start, xe, "conformal", cfg)
        num.append(tensor_power_from_mode(k, s.gamma[-1]))
        ex.append(power_exact(k, -xe / k))
    num, ex = np.array(num), np.array(ex)
    return dict(x_end=np.asarray(x_ends), numeric=num, exact=ex, rel_err=np.abs(num - ex) / ex,
                superhorizon=np.full_like(num, power_superhorizon()))


def sweep_H_Mpl(Hs=(0.5, 1.0, 2.0, 4.0), Ms=(1.0, 2.0, 4.0), k_over_H=1.0, x_start=100.0, x_end=1e-3):
    """Check the scaling Delta_t^2 = 2 H^2/(pi^2 M^2): wavenumber scaled with H."""
    rows = {"H": [], "M_pl": [], "numeric": [], "predicted": []}
    for H in Hs:
        for M in Ms:
            cfg = SolverConfig(H=H, M_pl=M, rtol=1e-11, atol=1e-14, n_out=10)
            k = k_over_H * H
            s = solve_mode(k, x_start, x_end, "conformal", cfg)
            rows["H"].append(H)
            rows["M_pl"].append(M)
            rows["numeric"].append(tensor_power_from_mode(k, s.gamma[-1]))
            rows["predicted"].append(power_exact(k, -x_end / k, H, M))
    return {k: np.array(v) for k, v in rows.items()}


def sweep_bogoliubov(phases=np.linspace(0, 2 * np.pi, 13), abs_beta=0.2, k=1.0, x_start=60.0, x_end=1e-4):
    """Superhorizon power for Bogoliubov-mixed initial states vs the prediction |alpha - beta|^2."""
    cfg = SolverConfig(rtol=1e-11, atol=1e-14, n_out=10)
    alpha = np.sqrt(1 + abs_beta**2)
    num, pred = [], []
    for ph in phases:
        beta = abs_beta * np.exp(1j * ph)
        s = solve_mode(k, x_start, x_end, "conformal", cfg, ic_name="bogoliubov", bogoliubov=(alpha, beta))
        ref = an.gamma_bd(k, s.eta[-1])
        num.append(abs(s.gamma[-1]) ** 2 / abs(ref) ** 2)
        u = ref
        pred.append(abs(alpha * u + beta * np.conj(u)) ** 2 / abs(u) ** 2)
    return dict(phase=np.asarray(phases), numeric=np.array(num), predicted=np.array(pred),
                superhorizon_leading=np.array([abs(alpha - abs_beta * np.exp(1j * p)) ** 2 for p in phases]))


def sweep_x_start(**kw):
    rows = cv.start_x_study(**kw)
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def sweep_tolerance(**kw):
    rows = cv.tolerance_study(**kw)
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}
