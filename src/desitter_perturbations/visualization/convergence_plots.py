"""Convergence figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S
from ..dynamics import convergence as cv


def plot_tolerance_and_trajectory():
    rows = cv.tolerance_study()
    traj = cv.error_trajectories()
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), constrained_layout=True)
    rt = np.array([r["rtol"] for r in rows]); me = np.array([r["max_err"] for r in rows])
    ax[0].loglog(rt, me, "o-", color=S.C["blue"], label="max relative error")
    ax[0].loglog(rt, rt, color=S.INK2, ls="--", lw=1.0, label="error = rtol")
    ax[0].set_xlabel("rtol (atol = 1e-3 rtol)"); ax[0].set_ylabel("max relative error vs exact")
    ax[0].set_title("DOP853: error tracks tolerance"); ax[0].legend(loc="upper left")
    for i, (rtol, (x, e)) in enumerate(traj.items()):
        ax[1].loglog(x[3:], e[3:] + 1e-18, color=S.SERIES[i], ls=S.LINESTYLES[i], label=f"rtol = {rtol:.0e}")
    ax[1].invert_xaxis(); ax[1].axvline(1.0, color=S.INK2, lw=1.0)
    ax[1].set_ylim(1e-11, 1e-2); S.direct_label(ax[1], 1.0, 3e-11, "x = 1", dx=-4, ha="right")
    ax[1].set_xlabel(r"$x=k/(aH)$"); ax[1].set_ylabel("relative error")
    ax[1].set_title("Error accumulates, no spike at horizon crossing"); ax[1].legend(loc="lower left", fontsize=8)
    return S.savefig(fig, "validation", "tolerance_and_horizon_crossing")


def plot_methods_and_order():
    meth = cv.method_comparison()
    h, e, p = cv.rk4_order_study()
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), constrained_layout=True)
    markers = ["o", "s", "^", "D", "v", "P"]
    for i, r in enumerate(meth):
        ax[0].loglog(r["nfev"], r["max_err"], markers[i], color=S.SERIES[i % 4] if i < 4 else S.C["magenta"], ms=7,
                     mec=S.SURFACE, label=r["method"])
    ax[0].set_xlabel("right-hand-side evaluations"); ax[0].set_ylabel("max relative error")
    ax[0].set_xticks([3e3, 1e4, 3e4, 1e5]); ax[0].set_xticklabels(["3e3", "1e4", "3e4", "1e5"]); ax[0].minorticks_off(); ax[0].set_title("Solver comparison at rtol = 1e-8"); ax[0].legend(loc="upper right", fontsize=8, ncol=2)
    ax[1].loglog(h, e, "o-", color=S.C["blue"], label=f"fixed-step RK4 (observed order {p:.2f})")
    ax[1].loglog(h, e[0] * (h / h[0]) ** 4, color=S.INK2, ls="--", lw=1.0, label=r"$h^4$ guide")
    ax[1].set_xticks([2e-3, 5e-3, 1e-2, 3e-2]); ax[1].set_xticklabels(["2e-3", "5e-3", "1e-2", "3e-2"]); ax[1].minorticks_off(); ax[1].set_xlabel(r"step in $N=\ln a$"); ax[1].set_ylabel("relative error at x = 0.01")
    ax[1].set_title("Resolution study (e-fold variable)"); ax[1].legend(loc="lower right", fontsize=8)
    return S.savefig(fig, "validation", "solver_methods_and_order"), p
