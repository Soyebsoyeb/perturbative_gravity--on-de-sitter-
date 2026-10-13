"""Levi-Civita connection, its first-order variation, and covariant derivatives
(component level, any coordinates)."""

from __future__ import annotations

import itertools
import sympy as sp


def christoffel(g: sp.Matrix, ginv: sp.Matrix, coords, simp=None):
    """Gamma[l][m][n] = (1/2) g^{lr} (d_m g_{rn} + d_n g_{rm} - d_r g_{mn})."""
    n = len(coords)
    dg = [[[sp.diff(g[i, j], coords[k]) for k in range(n)] for j in range(n)] for i in range(n)]
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for nn in range(m, n):
                s = 0
                for r in range(n):
                    if ginv[l, r] == 0:
                        continue
                    s += ginv[l, r] * (dg[r][nn][m] + dg[r][m][nn] - dg[m][nn][r])
                s = s / 2
                if simp:
                    s = simp(s)
                Gam[l][m][nn] = s
                Gam[l][nn][m] = s
    return Gam


def cov_deriv_covariant(T, Gam, coords, rank: int, simp=None):
    """Covariant derivative of a (0,rank) tensor stored as nested lists.

    Returns D with D[(a,) + idx] = nabla_a T_{idx}.
    """
    n = len(coords)
    D = {}
    for a in range(n):
        for idx in itertools.product(range(n), repeat=rank):
            val = sp.diff(_get(T, idx), coords[a])
            for pos in range(rank):
                for l in range(n):
                    new = list(idx)
                    new[pos] = l
                    val -= Gam[l][a][idx[pos]] * _get(T, tuple(new))
            if simp:
                val = simp(val)
            D[(a,) + idx] = val
    return D


def _get(T, idx):
    if isinstance(T, dict):
        return T[idx]
    if isinstance(T, sp.MatrixBase):
        return T[idx[0], idx[1]]
    out = T
    for i in idx:
        out = out[i]
    return out


def delta_christoffel(h: sp.Matrix, gbar_inv: sp.Matrix, Gam_bar, coords, simp=None):
    """First-order variation of the connection (a true tensor):

        dGamma^l_{mn} = (1/2) gbar^{lr} (nabla_m h_{rn} + nabla_n h_{rm} - nabla_r h_{mn}).
    """
    n = len(coords)
    Dh = cov_deriv_covariant(h, Gam_bar, coords, 2, simp=None)
    dG = [[[0] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for nn in range(n):
                s = 0
                for r in range(n):
                    if gbar_inv[l, r] == 0:
                        continue
                    s += gbar_inv[l, r] * (Dh[(m, r, nn)] + Dh[(nn, r, m)] - Dh[(r, m, nn)])
                s = s / 2
                dG[l][m][nn] = simp(s) if simp else s
    return dG


def delta_christoffel_direct(g_exact: sp.Matrix, ginv_exact: sp.Matrix, coords, eps, simp=None):
    """Brute-force dGamma = d/d eps Gamma[g_exact] at eps=0, with no covariant derivatives.

    Used as an independent check of :func:`delta_christoffel`.
    """
    Gam = christoffel(g_exact, ginv_exact, coords)
    n = len(coords)
    out = [[[0] * n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for nn in range(n):
                v = sp.diff(Gam[l][m][nn], eps).subs(eps, 0)
                out[l][m][nn] = simp(v) if simp else v
    return out
