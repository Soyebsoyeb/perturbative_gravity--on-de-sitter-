"""Einstein tensor, field equation E = G + Lambda g, and its linearization.

Two independent derivations of the first-order field equation are provided:

1. ``linearized_einstein_covariant``: background-covariant formulas
       dGamma^l_{mn} = 1/2 gbar^{lr}(nabla_m h_{rn} + nabla_n h_{rm} - nabla_r h_{mn})
       dR_{mn}       = nabla_l dGamma^l_{mn} - nabla_n dGamma^l_{l m}
       dR            = -h^{mn} Rbar_{mn} + gbar^{mn} dR_{mn}
       dG_{mn}       = dR_{mn} - 1/2 h_{mn} Rbar - 1/2 gbar_{mn} dR
       dE_{mn}       = dG_{mn} + Lambda h_{mn}
2. ``linearized_einstein_bruteforce``: component-level; builds g = gbar + eps h with a
   truncated inverse, computes G + Lambda g exactly to O(eps) and extracts the
   coefficient of eps.  No covariant derivative of h is used.
"""

from __future__ import annotations

import itertools
import sympy as sp

from .curvature import Geometry
from .connections import christoffel, delta_christoffel, cov_deriv_covariant
from .metric_expansion import inverse_metric_series, truncate, EPS


def einstein_tensor(geom: Geometry) -> sp.Matrix:
    return (geom.ricci - geom.g * geom.ricci_scalar / 2).applyfunc(geom.simp)


def field_equation(geom: Geometry, Lam) -> sp.Matrix:
    """E_{mn} = G_{mn} + Lambda g_{mn}."""
    return (einstein_tensor(geom) + Lam * geom.g).applyfunc(geom.simp)


def _cov_deriv_mixed_12(dG, Gam, coords, simp=None):
    """nabla_a T^l_{mn} for a (1,2) tensor T given as nested lists dG[l][m][n]."""
    n = len(coords)
    D = {}
    for a, l, m, nn in itertools.product(range(n), repeat=4):
        val = sp.diff(dG[l][m][nn], coords[a])
        for s in range(n):
            val += Gam[l][a][s] * dG[s][m][nn]
            val -= Gam[s][a][m] * dG[l][s][nn]
            val -= Gam[s][a][nn] * dG[l][m][s]
        D[(a, l, m, nn)] = simp(val) if simp else val
    return D


def linearized_einstein_covariant(geom_bar: Geometry, h: sp.Matrix, Lam, simp=None):
    """Return dict with dGamma, dRicci, dR, dG, dE for the covariant perturbation h."""
    simp = simp or geom_bar.simp
    n, coords = geom_bar.n, geom_bar.coords
    Gam = geom_bar.Gamma
    ginv = geom_bar.ginv
    dGam = delta_christoffel(h, ginv, Gam, coords, simp=simp)
    D = _cov_deriv_mixed_12(dGam, Gam, coords, simp=simp)
    dRic = sp.zeros(n, n)
    for m in range(n):
        for nn in range(n):
            # nabla_l dGamma^l_{mn} - nabla_n dGamma^l_{lm}
            v = sum(D[(l, l, m, nn)] for l in range(n)) - sum(D[(nn, l, l, m)] for l in range(n))
            dRic[m, nn] = simp(v)
    Ricbar = geom_bar.ricci
    Rbar = geom_bar.ricci_scalar
    h_up = ginv * h * ginv
    dR = simp(-sum(h_up[i, j] * Ricbar[i, j] for i in range(n) for j in range(n))
              + sum(ginv[i, j] * dRic[i, j] for i in range(n) for j in range(n)))
    dG = (dRic - h * Rbar / 2 - geom_bar.g * dR / 2).applyfunc(simp)
    dE = (dG + Lam * h).applyfunc(simp)
    return {"dGamma": dGam, "dRicci": dRic, "dR": dR, "dG": dG, "dE": dE}


def linearized_einstein_bruteforce(coords, gbar: sp.Matrix, h: sp.Matrix, Lam, simp=None):
    """Component-level dE_{mn} = coefficient of eps in G + Lam g for g = gbar + eps h."""
    n = len(coords)
    ginv_t = inverse_metric_series(gbar, h, order=1)
    g = (gbar + EPS * h)
    tr = lambda e: truncate(sp.expand(e), 1)
    # Christoffel truncated at O(eps)
    Gam = christoffel(g, ginv_t, coords)
    Gam = [[[tr(Gam[l][m][nn]) for nn in range(n)] for m in range(n)] for l in range(n)]
    Ric = sp.zeros(n, n)
    for s in range(n):
        for nn in range(s, n):
            val = 0
            for r in range(n):
                val += sp.diff(Gam[r][nn][s], coords[r]) - sp.diff(Gam[r][r][s], coords[nn])
                for l in range(n):
                    val += Gam[r][r][l] * Gam[l][nn][s] - Gam[r][nn][l] * Gam[l][r][s]
            val = tr(val)
            Ric[s, nn] = val
            Ric[nn, s] = val
    R = tr(sum(ginv_t[i, j] * Ric[i, j] for i in range(n) for j in range(n)))
    E = (Ric - g * R / 2 + Lam * g).applyfunc(tr)
    out0 = E.applyfunc(lambda e: e.coeff(EPS, 0))
    out1 = E.applyfunc(lambda e: e.coeff(EPS, 1))
    if simp:
        out0 = out0.applyfunc(simp)
        out1 = out1.applyfunc(simp)
    return out0, out1
