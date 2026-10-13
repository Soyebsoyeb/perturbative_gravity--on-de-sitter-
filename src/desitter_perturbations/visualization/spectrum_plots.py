"""Spectrum figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S
from ..spectra import parameter_sweeps as ps
from ..spectra.tensor_power import power_superhorizon, power_exact
from ..dynamics import convergence as cv


def plot_tensor_spectrum():
    r = ps.sweep_wavenumber(n=21, k_min=0.2, k_max=200.0, eta_end=-1e-2, x_start=300.0)
    kk = np.geomspace(0.1, 300, 400)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), constrained_layout=True)
    ax[0].loglog(kk, power_exact(kk, -1e-2), color=S.C["orange"], ls="--", lw=1.4, label=r"exact $\propto 1+k^2\eta^2$")
    ax[0].loglog(r["k"], r["numeric"], "o", color=S.C["blue"], ms=5, label="numerical modes")
    ax[0].axhline(power_superhorizon(), color=S.INK2, ls=":", lw=1.2)
    S.direct_label(ax[0], 0.12, power_superhorizon(), r"$2H^2/(\pi^2 M_{pl}^2)$", dx=0, dy=9)
    ax[0].axvline(100.0, color=S.C["aqua"], lw=1.0)
    S.direct_label(ax[0], 100.0, 2.0, r"$k|\eta_{end}|=1$", dx=-4, dy=0, ha="right")
    ax[0].set_xlabel("k / H"); ax[0].set_ylabel(r"$\Delta_t^2(k)$ at $\eta_{end}=-10^{-2}/H$")
    ax[0].set_title("Scale invariance outside the horizon, growth of still-subhorizon modes")
    ax[0].legend(loc="upper left", fontsize=8)
    ax[1].loglog(r["k"], r["rel_err"] + 1e-18, "o-", color=S.C["blue"], ms=4)
    ax[1].set_xlabel("k / H"); ax[1].set_ylabel("relative error vs exact")
    ax[1].set_title("Numerical vs analytic spectrum")
    return S.savefig(fig, "spectra", "tensor_spectrum")


def plot_late_time(k=1.0):
    r = ps.sweep_late_time()
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), constrained_layout=True)
    ax[0].loglog(r["x_end"], r["numeric"], "o", color=S.C["blue"], ms=4, label="numerical")
    ax[0].loglog(r["x_end"], r["exact"], color=S.C["orange"], ls="--", label=r"exact $(1+x^2)$")
    ax[0].axhline(power_superhorizon(), color=S.INK2, ls=":")
    ax[0].invert_xaxis()
    ax[0].set_xlabel(r"$x_{end}=k|\eta_{end}|$"); ax[0].set_ylabel(r"$\Delta_t^2$")
    ax[0].set_title("Freeze-out of the tensor spectrum"); ax[0].legend(loc="upper left")
    dev = r["numeric"] / power_superhorizon() - 1.0
    ax[1].loglog(r["x_end"], np.abs(dev), "o", color=S.C["blue"], ms=4, label=r"|$\Delta_t^2/\Delta_\infty^2 - 1$|")
    ax[1].loglog(r["x_end"], r["x_end"] ** 2, color=S.INK2, ls="--", lw=1.1, label=r"$x_{end}^2$")
    ax[1].invert_xaxis()
    ax[1].set_xlabel(r"$x_{end}$"); ax[1].set_title("Departure from scale invariance is $x^2$")
    ax[1].legend(loc="lower left")
    return S.savefig(fig, "spectra", "late_time_freezeout")


def plot_bogoliubov_and_start():
    b = ps.sweep_bogoliubov()
    rows = cv.start_x_study()
    xs = np.array([r["x_start"] for r in rows], float)
    dev = np.array([abs(r["ratio"] - 1) for r in rows])
    beta = np.array([r["abs_beta"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9), constrained_layout=True)
    ax[0].plot(b["phase"], b["numeric"], "o", color=S.C["blue"], ms=5, label="numerical")
    ax[0].plot(b["phase"], b["predicted"], color=S.C["orange"], ls="--", label="exact prediction")
    ax[0].plot(b["phase"], b["superhorizon_leading"], color=S.INK2, ls=":", label=r"leading $|\alpha-\beta|^2$")
    ax[0].set_xlabel(r"arg $\beta$"); ax[0].set_ylabel(r"$\Delta_t^2 / \Delta_{t,BD}^2$")
    ax[0].set_title(r"Excited initial state, $|\beta|=0.2$"); ax[0].legend(loc="lower center", fontsize=8)
    ax[1].loglog(xs, dev, "o", color=S.C["blue"], label="|P/P_BD - 1| (numerical)")
    ax[1].loglog(xs, 2 * beta, "s--", color=S.C["orange"], ms=4, label=r"$2|\beta|=1/x_{start}^2$")
    ax[1].set_xlabel(r"$x_{start}$"); ax[1].set_ylabel("error from leading-WKB initial data")
    ax[1].set_title("Finite-start-time error is a Bogoliubov mixing"); ax[1].legend(loc="lower left", fontsize=8)
    return S.savefig(fig, "spectra", "initial_state_dependence")


def plot_scaling_H_M():
    r = ps.sweep_H_Mpl()
    fig, ax = plt.subplots(figsize=(5.6, 4.2), constrained_layout=True)
    colors = {1.0: S.SERIES[0], 2.0: S.SERIES[1], 4.0: S.SERIES[2]}
    marks = {1.0: "o", 2.0: "s", 4.0: "^"}
    for M in (1.0, 2.0, 4.0):
        m = r["M_pl"] == M
        ax.loglog(r["H"][m], r["numeric"][m], marks[M], color=colors[M], label=f"M_pl = {M:g} (numerical)")
        ax.loglog(r["H"][m], r["predicted"][m], "-", color=colors[M], lw=1.0, alpha=0.7)
    ax.set_xlabel("H"); ax.set_ylabel(r"$\Delta_t^2$ at $k=H$, $x_{end}=10^{-3}$")
    ax.set_title(r"Scaling $\Delta_t^2 \propto H^2/M_{pl}^2$ (lines: prediction)"); ax.legend(loc="upper left", fontsize=8)
    return S.savefig(fig, "spectra", "scaling_H_Mpl")
