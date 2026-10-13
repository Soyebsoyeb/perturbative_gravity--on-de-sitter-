"""Gauge-comparison and SVT diagnostics figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S
from ..perturbations import gauge_invariants as gi
from ..perturbations import svt_numeric as sv


def plot_gauge_comparison():
    eta, res = gi.numeric_gauge_experiment(eta_range=(-50.0, -1.0))
    names = ["original", "shifted_1", "shifted_2"]
    labels = ["original gauge", "gauge shift 1", "gauge shift 2"]
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 6.6), constrained_layout=True)
    for i, n in enumerate(names):
        ax[0, 0].plot(eta, res[n]["A"], color=S.SERIES[i], ls=S.LINESTYLES[i], label=labels[i])
        ax[0, 1].plot(eta, res[n]["psi"], color=S.SERIES[i], ls=S.LINESTYLES[i], label=labels[i])
        ax[1, 0].plot(eta, res[n]["Phi"], color=S.SERIES[i], ls=S.LINESTYLES[i], label=labels[i])
    ax[0, 0].set_title("Gauge-dependent: lapse perturbation A")
    ax[0, 0].legend(loc="lower left")
    ax[0, 1].set_title(r"Gauge-dependent: curvature perturbation $\psi$")
    ax[1, 0].set_title(r"Invariant: Bardeen potential $\Phi_B$ (curves coincide)")
    for a in ax.flat[:3]:
        a.set_xlabel(r"$\eta$")
    # grouped bars: size of the change induced by each gauge shift
    quantities = [("A", r"$A$"), ("psi", r"$\psi$"), ("B", r"$B$"), ("E", r"$E$"), ("Phi", r"$\Phi_B$"), ("Psi", r"$\Psi_B$")]
    width = 0.38
    pos = np.arange(len(quantities))
    for j, shift in enumerate(("shifted_1", "shifted_2")):
        vals = [max(np.max(np.abs(res[shift][q] - res["original"][q])), 1e-18) for q, _ in quantities]
        ax[1, 1].bar(pos + (j - 0.5) * width, vals, width, color=S.SERIES[j + 1] if j == 0 else S.SERIES[2],
                     label=labels[j + 1])
    ax[1, 1].set_yscale("log")
    ax[1, 1].set_xticks(pos); ax[1, 1].set_xticklabels([l for _, l in quantities])
    ax[1, 1].set_ylim(1e-18, 1e2)
    ax[1, 1].set_ylabel("max |change| over eta")
    ax[1, 1].set_title("Change under the gauge shift: O(1) vs round-off")
    ax[1, 1].legend(loc="center right")
    return S.savefig(fig, "gauge", "gauge_comparison")


def plot_svt_diagnostics():
    Ns = [8, 12, 16, 20, 24]
    rec = [sv.reconstruction_error(sv.random_symmetric_field(N, seed=2)) for N in Ns]
    spec, fd = [], []
    for N in Ns:
        g = sv.tt_part_field(N, seed=4)
        spec.append(np.max(np.abs(sv.spectral_divergence(g))) + 1e-20)
        fd.append(np.max(np.abs(sv.finite_difference_divergence(g))))
    kn, eE, eF = sv.noise_amplification()
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8), constrained_layout=True)
    ax[0].semilogy(Ns, np.maximum(rec, 1e-19), "o-", color=S.C["blue"])
    ax[0].set_xlabel("grid points per side N"); ax[0].set_ylabel("max |reconstruction - h|")
    ax[0].set_title("SVT round trip (spectral)"); ax[0].set_ylim(1e-19, 1e-14)
    ax[1].semilogy(Ns, spec, "o-", color=S.C["aqua"], label="spectral divergence")
    ax[1].semilogy(Ns, fd, "s--", color=S.C["orange"], label="2nd-order finite difference")
    ax[1].set_xlabel("N"); ax[1].set_ylabel(r"max $|\partial_j\gamma_{ij}|$")
    ax[1].set_title("Transversality of the TT part"); ax[1].legend(loc="center right")
    ax[2].loglog(kn, eE, ".", color=S.C["blue"], ms=3, alpha=0.5, label="error in E")
    ax[2].loglog(kn, eF, ".", color=S.C["orange"], ms=3, alpha=0.5, label="error in F")
    kk = np.linspace(kn.min(), kn.max(), 50)
    ax[2].loglog(kk, 1e-8 * kk**-2 * 1.0, color=S.INK2, lw=1, ls="--")
    ax[2].loglog(kk, 1e-8 * kk**-1 * 1.0, color=S.INK2, lw=1, ls=":")
    ax[2].set_xlabel("|k|"); ax[2].set_ylabel("recovery error for noise 1e-8")
    ax[2].set_title("Small-k amplification: guides $k^{-2}$, $k^{-1}$"); ax[2].legend(loc="lower left")
    return S.savefig(fig, "perturbations", "svt_diagnostics")
