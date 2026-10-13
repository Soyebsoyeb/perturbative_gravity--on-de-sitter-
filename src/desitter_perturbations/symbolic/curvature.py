"""Riemann, Ricci and scalar curvature from a metric given in explicit coordinates."""

from __future__ import annotations

import itertools
import sympy as sp

from .connections import christoffel
from ..conventions import RIEMANN_SIGN


class Geometry:
    """Curvature of a metric g(coords).

    Parameters
    ----------
    coords : sequence of sympy symbols
    g : sympy Matrix (covariant metric)
    ginv : optional inverse (for truncated perturbative inverses)
    simp : optional callable applied to every component (default: sp.simplify)
    sign : +1 is the convention of docs/conventions.md; -1 flips the Riemann
           sign and exists only so that tests can demonstrate sign detection.
    """

    def __init__(self, coords, g, ginv=None, simp=sp.simplify, sign: int = RIEMANN_SIGN):
        self.coords = tuple(coords)
        self.n = len(self.coords)
        self.g = sp.Matrix(g)
        self.ginv = sp.Matrix(ginv) if ginv is not None else self.g.inv()
        self.simp = simp if simp is not None else (lambda e: e)
        self.sign = sign
        self._Gam = None
        self._Riem_up = None
        self._Ric = None
        self._R = None

    # ---------------------------------------------------------------- pieces
    @property
    def Gamma(self):
        if self._Gam is None:
            self._Gam = christoffel(self.g, self.ginv, self.coords, simp=self.simp)
        return self._Gam

    @property
    def riemann_up(self):
        """R^r_{s m n} as a dict keyed by (r,s,m,n)."""
        if self._Riem_up is None:
            n, c, G = self.n, self.coords, self.Gamma
            R = {}
            for r, s, m, nn in itertools.product(range(n), repeat=4):
                if m == nn:
                    R[(r, s, m, nn)] = 0
                    continue
                if m > nn:
                    R[(r, s, m, nn)] = -R[(r, s, nn, m)]
                    continue
                val = sp.diff(G[r][nn][s], c[m]) - sp.diff(G[r][m][s], c[nn])
                for l in range(n):
                    val += G[r][m][l] * G[l][nn][s] - G[r][nn][l] * G[l][m][s]
                R[(r, s, m, nn)] = self.simp(self.sign * val)
            self._Riem_up = R
        return self._Riem_up

    def riemann_down(self):
        """R_{r s m n}."""
        n = self.n
        Ru = self.riemann_up
        out = {}
        for r, s, m, nn in itertools.product(range(n), repeat=4):
            out[(r, s, m, nn)] = self.simp(sum(self.g[r, l] * Ru[(l, s, m, nn)] for l in range(n)))
        return out

    @property
    def ricci(self) -> sp.Matrix:
        """R_{s n} = R^r_{s r n} (uses the direct contraction, cheaper than full Riemann)."""
        if self._Ric is None:
            n, c, G = self.n, self.coords, self.Gamma
            Ric = sp.zeros(n, n)
            for s in range(n):
                for nn in range(s, n):
                    val = 0
                    for r in range(n):
                        val += sp.diff(G[r][nn][s], c[r]) - sp.diff(G[r][r][s], c[nn])
                        for l in range(n):
                            val += G[r][r][l] * G[l][nn][s] - G[r][nn][l] * G[l][r][s]
                    val = self.simp(self.sign * val)
                    Ric[s, nn] = val
                    Ric[nn, s] = val
            self._Ric = Ric
        return self._Ric

    @property
    def ricci_scalar(self):
        if self._R is None:
            self._R = self.simp(sum(self.ginv[i, j] * self.ricci[i, j]
                                    for i in range(self.n) for j in range(self.n)))
        return self._R

    # ------------------------------------------------------------ invariants
    def kretschmann(self):
        """R_{abcd} R^{abcd}."""
        n = self.n
        Rd = self.riemann_down()
        # raise all indices with ginv
        tot = 0
        for a, b, c_, d in itertools.product(range(n), repeat=4):
            if Rd[(a, b, c_, d)] == 0:
                continue
            up = 0
            for a2, b2, c2, d2 in itertools.product(range(n), repeat=4):
                if self.ginv[a, a2] == 0 or self.ginv[b, b2] == 0 or self.ginv[c_, c2] == 0 or self.ginv[d, d2] == 0:
                    continue
                up += (self.ginv[a, a2] * self.ginv[b, b2] * self.ginv[c_, c2] *
                       self.ginv[d, d2] * Rd[(a2, b2, c2, d2)])
            tot += Rd[(a, b, c_, d)] * up
        return self.simp(tot)

    def ricci_squared(self):
        Ric_up = self.ginv * self.ricci * self.ginv
        return self.simp(sum(self.ricci[i, j] * Ric_up[i, j]
                             for i in range(self.n) for j in range(self.n)))

    def constant_curvature_residual(self, H):
        """Max-simplified residual of R_{mnrs} - H^2 (g_{mr} g_{ns} - g_{ms} g_{nr})."""
        n = self.n
        Rd = self.riemann_down()
        res = {}
        for a, b, c_, d in itertools.product(range(n), repeat=4):
            target = H**2 * (self.g[a, c_] * self.g[b, d] - self.g[a, d] * self.g[b, c_])
            res[(a, b, c_, d)] = sp.simplify(Rd[(a, b, c_, d)] - target)
        return res
