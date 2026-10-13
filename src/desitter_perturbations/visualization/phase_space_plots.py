"""Phase-space and conserved-quantity figures."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from . import style as S
from ..dynamics import analytic_solutions as an
from ..dynamics.numerical_solvers import SolverConfig, solve_mode
from ..quantization import wronskian as wr


def plot_phase_space(k=1.0):
    cfg = SolverConfig(rtol=1e-11, atol=1e-14, n_out=2500)
    s = solve_mode(k, 60.0, 1e-2, "conformal", cfg)
    z = an.z_T(s.eta)
    v = z * s.gamma
    dv = z * (s.dgamma_deta + (-1.0 / s.eta) * s.gamma)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.0), constrained_layout=True)
    sc = ax[0].scatter(np.real(v) * np.sqrt(2 * k), np.real(dv) * np.sqrt(2 * k) / k, c=np.log10(s.x),
                       s=3, cmap="Blues_r", vmin=-2, vmax=2)
    ax[0].set_xlabel(r"Re $v\,\sqrt{2k}$"); ax[0].set_ylabel(r"Re $v'\,\sqrt{2k}/k$")
    ax[0].set_title("Canonical phase space (colour: log10 x)")
    ax[0].set_xlim(-2.2, 2.2); ax[0].set_ylim(-2.2, 2.2); ax[0].set_aspect('equal')
    fig.colorbar(sc, ax=ax[0], label="log10 x", shrink=0.8)
    ax[1].semilogx(s.x, np.real(v) * np.sqrt(2 * k), color=S.C["blue"], label="Re v")
    ax[1].semilogx(s.x, np.imag(v) * np.sqrt(2 * k), color=S.C["orange"], ls="--", label="Im v")
    ax[1].invert_xaxis(); ax[1].set_xlabel("x"); ax[1].set_title(r"Mode $v$: oscillation, then growth")
    ax[1].set_ylim(-12, 12); ax[1].legend(loc="upper left")
    W = np.abs(wr.gamma_wronskian(s.gamma, s.dgamma_deta, z) - 1j)
    ax[2].loglog(s.x[3:], W[3:] + 1e-18, color=S.C["blue"])
    ax[2].invert_xaxis(); ax[2].set_xlabel("x"); ax[2].set_ylabel(r"$|W - i|$")
    ax[2].set_title("Wronskian conservation (commutator)")
    return S.savefig(fig, "dynamics", "phase_space_and_wronskian")
