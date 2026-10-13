"""Background geometry figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S
from .. import background as bg


def plot_background_evolution(H=1.0):
    t = np.linspace(-2.0, 6.0, 600) / H
    a = bg.scale_factor_t(t, H)
    eta = bg.conformal_time_from_t(t, H)
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.4), constrained_layout=True)
    ax[0].semilogy(t * H, a, color=S.C["blue"])
    ax[0].set_xlabel("H t"); ax[0].set_ylabel("a(t)"); ax[0].set_title("Scale factor a = exp(Ht)")
    ax[1].plot(t * H, -eta * H, color=S.C["orange"])
    ax[1].set_yscale("log")
    ax[1].set_xlabel("H t"); ax[1].set_ylabel(r"$-H\eta$"); ax[1].set_title(r"Conformal time $\eta=-e^{-Ht}/H$")
    hub = bg.hubble_radius_comoving(eta)
    ax[2].semilogy(t * H, hub * H, color=S.C["aqua"], label=r"comoving Hubble radius $1/(aH)=-\eta$")
    ax[2].semilogy(t * H, 1.0 / (a * H) * H, color=S.INK2, ls=":", lw=1.2, label=r"$1/(aH)$ (same curve)")
    ax[2].set_xlabel("H t"); ax[2].set_ylabel(r"$H/(aH)\ [1/H]$"); ax[2].set_title("Shrinking comoving horizon")
    ax[2].legend(loc="upper right")
    return S.savefig(fig, "geometry", "background_evolution")


def plot_horizon_crossing(H=1.0, ks=(1.0, 10.0, 100.0)):
    N = np.linspace(-1.0, 9.0, 800)           # e-folds (a = e^N)
    a = np.exp(N)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9), constrained_layout=True)
    for i, k in enumerate(ks):
        ax[0].semilogy(N, k / a / H, color=S.SERIES[i], ls=S.LINESTYLES[i], label=f"k/H = {k:g}")
        Nc = np.log(k / H)
        ax[0].plot([Nc], [1.0], "o", color=S.SERIES[i], mec=S.SURFACE, mew=1.5)
    ax[0].axhline(1.0, color=S.INK2, lw=1.0)
    ax[0].set_xlabel("N = ln a"); ax[0].set_ylabel(r"physical wavenumber $k/(aH)$")
    ax[0].set_title("Horizon crossing at k/(aH) = 1")
    ax[0].text(8.9, 1.25, "Hubble scale", ha="right", color=S.INK2, fontsize=8.5)
    ax[0].legend(loc="upper right", bbox_to_anchor=(1.0, 0.9))
    # conformal diagram of comoving scales
    eta = -np.geomspace(1e-3, 10, 600)
    for i, k in enumerate(ks):
        ax[1].loglog(-eta, np.full_like(eta, 1.0 / k), color=S.SERIES[i], ls=S.LINESTYLES[i])
    ax[1].loglog(-eta, -eta, color=S.INK, lw=1.4)
    ax[1].invert_xaxis()
    ax[1].set_xlabel(r"$-\eta H$ (earlier at left)"); ax[1].set_ylabel(r"comoving length $\times H$")
    ax[1].set_title(r"Comoving wavelengths $1/k$ vs Hubble radius $-\eta$")
    S.direct_label(ax[1], 0.0016, 0.0016, "Hubble radius", dx=6, dy=10)
    return S.savefig(fig, "geometry", "horizon_crossing")
