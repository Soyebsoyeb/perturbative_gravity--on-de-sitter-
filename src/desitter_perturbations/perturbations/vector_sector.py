"""Vector sector: h_0x = a^2 S_x(eta,z), h_xz = a^2 d_z F_x(eta,z).

Gauge invariant V = S - F'.  Verified exactly:

    dE_0x = -(1/2) d_z^2 V,        dE_xz = -(1/2) d_z (V' + 2 calH V).
"""

from __future__ import annotations

import functools

import sympy as sp

from ._common import setup, linearized_dE


@functools.lru_cache(maxsize=None)
def vector_dE():
    s = setup()
    eta, z, a, cH = s["eta"], s["z"], s["a"], s["calH"]
    S = sp.Function("S")(eta, z)
    F = sp.Function("F")(eta, z)
    h = sp.zeros(4, 4)
    h[0, 1] = h[1, 0] = a**2 * S
    h[1, 3] = h[3, 1] = a**2 * sp.diff(F, z)
    dE = linearized_dE(s, h)
    V = S - sp.diff(F, eta)
    res = {
        "0x": sp.simplify(dE[0, 1] + sp.diff(V, z, 2) / 2),
        "xz": sp.simplify(dE[1, 3] + sp.diff(sp.diff(V, eta) + 2 * cH * V, z) / 2),
    }
    return {"dE": dE, "V": V, "residuals": res, "setup": s}
