"""Linearized Bianchi identity.

The exact identity nabla_m G^m_n = 0 holds for every metric.  Linearizing around
a background with nabla-bar Gbar^m_n = 0 gives (mixed-index form)

    nabla-bar_m dG^m_n + dGamma^m_{m l} Gbar^l_n - dGamma^l_{m n} Gbar^m_l = 0,
    dG^m_n = gbar^{m r} dG_{r n} - h^{m r} Gbar_{r n}.

For the covariant-index form the same statement reads

    gbar^{m r} [ nabla-bar_r dG_{m n} - dGamma^l_{r m} Gbar_{l n} - dGamma^l_{r n} Gbar_{m l} ] = 0,

i.e. the naive divergence nabla-bar^m dG_{mn} does NOT vanish by itself off the
background field equations: the variation of the connection contributes terms
proportional to Gbar = -Lambda gbar.  For the field equation E = G + Lambda g
(which vanishes on the background) these extra terms are absent and

    nabla-bar^m dE_{mn} = 0

holds exactly.  All three statements are checked componentwise.
"""

from __future__ import annotations

import itertools
import sympy as sp

from .connections import cov_deriv_covariant


def _div_covariant_2(T: sp.Matrix, Gam, ginv, coords, simp):
    """gbar^{m r} nabla_r T_{m n}."""
    n = len(coords)
    D = cov_deriv_covariant(T, Gam, coords, 2)
    return [simp(sum(ginv[m, r] * D[(r, m, nu)] for m in range(n) for r in range(n))) for nu in range(n)]


def bianchi_covariant_residual(geom_bar, dG: sp.Matrix, dGam, simp=None):
    """Residual of the covariant-index linearized Bianchi identity for dG."""
    simp = simp or geom_bar.simp
    n = geom_bar.n
    Gam, ginv = geom_bar.Gamma, geom_bar.ginv
    Gbar = (geom_bar.ricci - geom_bar.g * geom_bar.ricci_scalar / 2)
    naive = _div_covariant_2(dG, Gam, ginv, geom_bar.coords, simp)
    corr = []
    for nu in range(n):
        s = 0
        for m, r, l in itertools.product(range(n), repeat=3):
            s += ginv[m, r] * (dGam[l][r][m] * Gbar[l, nu] + dGam[l][r][nu] * Gbar[m, l])
        corr.append(simp(s))
    return [simp(naive[i] - corr[i]) for i in range(n)], naive, corr


def bianchi_mixed_residual(geom_bar, h: sp.Matrix, dG: sp.Matrix, dGam, simp=None):
    """Residual of the mixed-index linearized Bianchi identity."""
    simp = simp or geom_bar.simp
    n = geom_bar.n
    coords, Gam, ginv = geom_bar.coords, geom_bar.Gamma, geom_bar.ginv
    Gbar_dn = (geom_bar.ricci - geom_bar.g * geom_bar.ricci_scalar / 2)
    Gbar_mixed = ginv * Gbar_dn  # G^m_n
    dGmix = ginv * dG - (ginv * h * ginv) * Gbar_dn  # dG^m_n
    res = []
    for nu in range(n):
        s = 0
        for m in range(n):
            d = sp.diff(dGmix[m, nu], coords[m])
            for l in range(n):
                d += Gam[m][m][l] * dGmix[l, nu] - Gam[l][m][nu] * dGmix[m, l]
            s += d
        for m, l in itertools.product(range(n), repeat=2):
            s += dGam[m][m][l] * Gbar_mixed[l, nu] - dGam[l][m][nu] * Gbar_mixed[m, l]
        res.append(simp(s))
    return res


def bianchi_field_equation_residual(geom_bar, dE: sp.Matrix, simp=None):
    """nabla-bar^m dE_{mn}; must vanish identically on a background solution."""
    simp = simp or geom_bar.simp
    return _div_covariant_2(dE, geom_bar.Gamma, geom_bar.ginv, geom_bar.coords, simp)
