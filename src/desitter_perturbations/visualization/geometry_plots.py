"""Geometry and perturbative-expansion figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S


def _random_background_and_h(rng, scale=1.0):
    d = np.array([-1.3, 0.8, 1.1, 0.6])
    gbar = np.diag(d)
    A = rng.standard_normal((4, 4))
    h = scale * 0.5 * (A + A.T)
    return gbar, h


def expansion_errors(eps_values, seed=7):
    rng = np.random.default_rng(seed)
    gbar, h = _random_background_and_h(rng)
    gi = np.linalg.inv(gbar)
    M = gi @ h
    out = {"inv0": [], "inv1": [], "inv2": [], "det0": [], "det1": [], "det2": []}
    tr = np.trace(M)
    tr2 = np.trace(M @ M)
    for e in eps_values:
        g = gbar + e * h
        ginv_exact = np.linalg.inv(g)
        t0 = gi
        t1 = gi - e * M @ gi
        t2 = t1 + e**2 * M @ M @ gi
        for name, t in (("inv0", t0), ("inv1", t1), ("inv2", t2)):
            out[name].append(np.max(np.abs(t - ginv_exact)))
        ratio = np.sqrt(abs(np.linalg.det(g) / np.linalg.det(gbar)))
        s0, s1, s2 = 1.0, 1 + e * tr / 2, 1 + e * tr / 2 + e**2 * (tr**2 / 8 - tr2 / 4)
        for name, s in (("det0", s0), ("det1", s1), ("det2", s2)):
            out[name].append(abs(s - ratio))
    return {k: np.array(v) for k, v in out.items()}


def plot_expansion_errors():
    eps = np.geomspace(1e-3, 3e-2, 14)
    r = expansion_errors(eps)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9), constrained_layout=True)
    slopes = {}
    for panel, (title, keys) in enumerate((("Inverse metric truncation error", ("inv0", "inv1", "inv2")),
                                           (r"$\sqrt{-g}$ truncation error", ("det0", "det1", "det2")))):
        for i, key in enumerate(keys):
            ax[panel].loglog(eps, r[key], marker="o", color=S.SERIES[i], ls=S.LINESTYLES[i], ms=4,
                             label=f"through O($\\epsilon^{i}$)")
            slopes[key] = float(np.polyfit(np.log(eps), np.log(r[key]), 1)[0])
        ax[panel].set_xlabel(r"$\epsilon$"); ax[panel].set_ylabel("max abs error")
        ax[panel].set_title(title)
        ax[panel].legend(loc="lower right")
    for i, key in enumerate(("inv0", "inv1", "inv2")):
        S.direct_label(ax[0], eps[-1], r[key][-1], f"slope {slopes[key]:.2f}", dx=5)
    for i, key in enumerate(("det0", "det1", "det2")):
        S.direct_label(ax[1], eps[-1], r[key][-1], f"slope {slopes[key]:.2f}", dx=5)
    return S.savefig(fig, "geometry", "expansion_truncation_errors"), slopes


def plot_curvature_invariants():
    """Ricci scalar from the symbolic engine in both slicings, for the documented
    Riemann sign (+1) and a deliberately flipped sign (-1), against the analytic 12 H^2.

    The flipped-sign bars show that the convention test has teeth: a sign error is
    visible immediately as R = -12 H^2.
    """
    from ..background import symbolic_background
    from ..symbolic.curvature import Geometry
    b = symbolic_background()
    H = b["H"]
    vals = {}
    for slicing, key, coords in (("cosmic t", "g_cosmic", b["coords_cosmic"]),
                                 ("conformal eta", "g_conformal", b["coords_conformal"])):
        for sign in (+1, -1):
            R = Geometry(coords, b[key], sign=sign).ricci_scalar
            vals[(slicing, sign)] = float(R / H**2)
    fig, ax = plt.subplots(figsize=(6.6, 3.7), constrained_layout=True)
    pos = np.arange(2)
    ax.bar(pos - 0.18, [vals[("cosmic t", 1)], vals[("conformal eta", 1)]], 0.34,
           color=S.C["blue"], label="documented sign (+1)")
    ax.bar(pos + 0.18, [vals[("cosmic t", -1)], vals[("conformal eta", -1)]], 0.34,
           color=S.C["orange"], label="flipped Riemann sign (-1): detected")
    ax.axhline(12, color=S.INK2, lw=1.0, ls="--")
    S.direct_label(ax, 1.45, 12, "analytic 12", dx=0, dy=7, ha="right")
    ax.set_xticks(pos); ax.set_xticklabels(["cosmic-time slicing", "conformal-time slicing"])
    ax.set_ylabel(r"$R/H^2$ from the engine"); ax.set_title("Curvature-sign convention check")
    ax.legend(loc="lower right")
    return S.savefig(fig, "geometry", "curvature_sign_check")
