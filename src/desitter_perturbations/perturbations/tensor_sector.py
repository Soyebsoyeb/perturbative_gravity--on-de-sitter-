"""Tensor sector: TT field equation, polarization tensors, TT projection.

Symbolic result (verified): for h_ij = a^2 gamma_ij with gamma_xx = -gamma_yy = f_p,
gamma_xy = f_c depending on (eta, z),

    dE_ij = (1/2) (gamma_ij'' + 2 calH gamma_ij' - nabla^2 gamma_ij),

with all other components zero.
"""

from __future__ import annotations

import functools

import numpy as np
import sympy as sp

from ._common import setup, linearized_dE


@functools.lru_cache(maxsize=None)
def tensor_dE():
    s = setup()
    eta, z, a, cH = s["eta"], s["z"], s["a"], s["calH"]
    fp = sp.Function("f_p")(eta, z)
    fc = sp.Function("f_c")(eta, z)
    gam = sp.Matrix([[0, 0, 0, 0], [0, fp, fc, 0], [0, fc, -fp, 0], [0, 0, 0, 0]])
    h = a**2 * gam
    dE = linearized_dE(s, h)
    box = lambda f: sp.diff(f, eta, 2) + 2 * cH * sp.diff(f, eta) - sp.diff(f, z, 2)
    res = {
        "xx": sp.simplify(dE[1, 1] - box(fp) / 2),
        "yy": sp.simplify(dE[2, 2] + box(fp) / 2),
        "xy": sp.simplify(dE[1, 2] - box(fc) / 2),
        "others": [sp.simplify(dE[i, j]) for (i, j) in [(0, 0), (0, 1), (0, 2), (0, 3), (1, 3), (2, 3), (3, 3)]],
    }
    return {"dE": dE, "residuals": res, "setup": s}


# ------------------------------------------------------------ numeric tools

def transverse_projector(k):
    k = np.asarray(k, dtype=float)
    k2 = k @ k
    if k2 == 0.0:
        raise ValueError("projector is undefined at k = 0")
    return np.eye(3) - np.outer(k, k) / k2


def tt_projector(k):
    """Lambda_{ij,kl}(k) = 1/2 (P_ik P_jl + P_il P_jk - P_ij P_kl)."""
    P = transverse_projector(k)
    return 0.5 * (np.einsum("ik,jl->ijkl", P, P) + np.einsum("il,jk->ijkl", P, P)
                  - np.einsum("ij,kl->ijkl", P, P))


def tt_project(h, k):
    return np.einsum("ijkl,kl->ij", tt_projector(k), h)


def polarization_basis(k):
    """Real linear polarization tensors e^+, e^x with e^lam_ij e^lam'_ij = 2 delta."""
    k = np.asarray(k, dtype=float)
    khat = k / np.linalg.norm(k)
    trial = np.array([1.0, 0.0, 0.0]) if abs(khat[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(khat, trial)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(khat, e1)
    e_plus = np.outer(e1, e1) - np.outer(e2, e2)
    e_cross = np.outer(e1, e2) + np.outer(e2, e1)
    return e_plus, e_cross


@functools.lru_cache(maxsize=None)
def covariant_tensor_equation():
    """Verify the background-covariant form of the TT equation on de Sitter space.

    For h_mn = a^2 gamma_mn with the plus/cross ansatz (transverse, traceless):

        nabla^m h_mn = 0,   h = 0,
        dE_mn = -(1/2) (box h_mn - 2 H^2 h_mn),    box = gbar^{ab} nabla_a nabla_b,

    i.e. the massless spin-2 equation on de Sitter space (Lichnerowicz operator form).
    Returns residual dict (all entries should be zero).
    """
    from ..symbolic.connections import cov_deriv_covariant
    s = setup()
    eta, z, a, H, geom = s["eta"], s["z"], s["a"], s["H"], s["geom"]
    fp = sp.Function("f_p")(eta, z)
    fc = sp.Function("f_c")(eta, z)
    gam = sp.Matrix([[0, 0, 0, 0], [0, fp, fc, 0], [0, fc, -fp, 0], [0, 0, 0, 0]])
    h = a**2 * gam
    Gam, ginv, coords = geom.Gamma, geom.ginv, geom.coords
    D1 = cov_deriv_covariant(h, Gam, coords, 2)
    div = [sp.simplify(sum(ginv[m, r] * D1[(r, m, n)] for m in range(4) for r in range(4))) for n in range(4)]
    tr = sp.simplify(sum(ginv[m, n] * h[m, n] for m in range(4) for n in range(4)))
    D2 = cov_deriv_covariant(D1, Gam, coords, 3)       # D2[(a, b, m, n)] = nabla_a nabla_b h_mn
    box = sp.Matrix(4, 4, lambda m, n: sp.simplify(sum(ginv[p, q] * D2[(p, q, m, n)]
                                                       for p in range(4) for q in range(4))))
    dE = linearized_dE(s, h)
    resid = (dE + (box - 2 * H**2 * h) / 2).applyfunc(sp.simplify)
    return {"divergence": div, "trace": tr, "equation": resid}
