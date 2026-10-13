"""Tensor-mode evolution figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S
from ..dynamics import analytic_solutions as an
from ..dynamics import convergence as cv
from ..dynamics import residuals as rs
from ..dynamics.numerical_solvers import SolverConfig, solve_mode


def plot_mode_evolution(k=1.0):
    cfg = SolverConfig(rtol=1e-10, atol=1e-13, n_out=1500)
    s = solve_mode(k, 200.0, 1e-3, "conformal", cfg)
    ref = an.gamma_ref(k)
    ex = an.gamma_bd(k, s.eta) / ref
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), constrained_layout=True)
    ax[0].loglog(s.x, np.abs(s.gamma) / ref, color=S.C["blue"], label="numerical |gamma|")
    ax[0].loglog(s.x, np.abs(ex), color=S.C["orange"], ls="--", lw=1.4, label=r"exact $\sqrt{1+x^2}$")
    ax[0].loglog(s.x, s.x, color=S.INK2, ls=":", lw=1.2, label="subhorizon asymptote x")
    ax[0].axhline(1.0, color=S.INK2, ls="-.", lw=1.0)
    ax[0].axvline(1.0, color=S.C["aqua"], lw=1.2)
    S.direct_label(ax[0], 1.0, 8e-1, "horizon crossing x = 1", dx=6, dy=-10)
    S.direct_label(ax[0], 3e-2, 1.0, "superhorizon: frozen", dx=0, dy=9, ha="center")
    ax[0].invert_xaxis()
    ax[0].set_xlabel(r"$x = k/(aH)$ (time runs to the right)")
    ax[0].set_ylabel(r"$|\gamma_k|\, M_{pl} k^{3/2}/H$")
    ax[0].set_title("Amplitude: decay inside, freeze-out outside the horizon")
    ax[0].legend(loc="upper center", fontsize=8)
    m = s.x < 30
    ax[1].semilogx(s.x[m], np.real(s.gamma[m] / ref), color=S.C["blue"], label="Re gamma (numerical)")
    ax[1].semilogx(s.x[m], np.real(ex[m]), color=S.C["orange"], ls="--", lw=1.2, label="Re gamma (exact)")
    ax[1].invert_xaxis()
    ax[1].set_xlabel(r"$x = k/(aH)$"); ax[1].set_ylabel(r"Re $\gamma_k$ / reference")
    ax[1].set_title("Oscillation stops at horizon exit"); ax[1].legend(loc="upper left", fontsize=8)
    return S.savefig(fig, "dynamics", "mode_evolution")


def plot_formulation_comparison(k=1.0):
    sols, errs = cv.formulation_comparison(k=k)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), constrained_layout=True)
    names = {"conformal": "conformal time", "cosmic": "cosmic time", "efold": "e-folds N",
             "canonical": "canonical v"}
    for i, (f, s) in enumerate(sols.items()):
        ax[0].loglog(s.x[3:], rs.error_vs_exact(s)[3:] + 1e-18, color=S.SERIES[i], ls=S.LINESTYLES[i], label=names[f])
    ax[0].invert_xaxis()
    ax[0].axvline(1.0, color=S.INK2, lw=1.0)
    ax[0].set_xlabel(r"$x = k/(aH)$"); ax[0].set_ylabel("relative error vs exact BD")
    ax[0].set_ylim(1e-11, 1e-8); ax[0].set_title("Four formulations, same tolerance"); ax[0].legend(loc="lower left", fontsize=8)
    ref = sols["conformal"]
    for i, f in enumerate(("cosmic", "efold", "canonical")):
        d = np.abs(sols[f].gamma - ref.gamma) / np.abs(ref.gamma) + 1e-18
        ax[1].loglog(ref.x[3:], d[3:], color=S.SERIES[i + 1], ls=S.LINESTYLES[i + 1], label=f"{names[f]} vs conformal")
    ax[1].invert_xaxis()
    ax[1].set_xlabel(r"$x = k/(aH)$"); ax[1].set_ylabel("pointwise relative difference")
    ax[1].set_ylim(1e-12, 1e-8); ax[1].set_title("Cross-formulation agreement"); ax[1].legend(loc="lower left", fontsize=8)
    return S.savefig(fig, "dynamics", "formulation_comparison")


def plot_residual_diagnostics(k=1.0):
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), constrained_layout=True)
    ns = [1000, 2000, 4000, 8000]
    peak = []
    for i, n in enumerate(ns):
        s = solve_mode(k, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-12, atol=1e-15, n_out=n))
        kin, dyn = rs.equation_residuals(s)
        mask = s.x < 50
        peak.append(float(np.max(dyn[mask])))
        ax[0].loglog(s.x[mask], dyn[mask] + 1e-18, color=S.SERIES[i], ls=S.LINESTYLES[i], label=f"{n} output points")
    ax[0].invert_xaxis()
    ax[0].set_xlabel(r"$x$"); ax[0].set_ylabel("normalised ODE residual")
    ax[0].set_title("Residual of the stored solution (spline derivative)"); ax[0].legend(loc="lower left", fontsize=8)
    ax[1].loglog(ns, peak, "o-", color=S.C["blue"])
    ref = np.array(ns, float)
    ax[1].loglog(ref, peak[0] * (ref / ref[0]) ** -4.0, color=S.INK2, ls="--", lw=1.1, label=r"slope $-4$ guide")
    ax[1].set_xlabel("output points"); ax[1].set_ylabel("peak residual")
    ax[1].set_title("Residual set by interpolation error, not by the solver"); ax[1].legend(loc="upper right")
    return S.savefig(fig, "validation", "residual_diagnostics")
