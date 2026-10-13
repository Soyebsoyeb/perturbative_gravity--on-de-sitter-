"""Shared symbolic setup for the sector modules (conformal-time flat de Sitter)."""

from __future__ import annotations

import sympy as sp

from ..background import symbolic_background
from ..symbolic.curvature import Geometry
from ..symbolic.einstein_tensor import linearized_einstein_covariant


def setup():
    b = symbolic_background()
    H, eta = b["H"], b["eta"]
    x, y, z = b["spatial"]
    geom = Geometry(b["coords_conformal"], b["g_conformal"])
    return {
        "b": b, "H": H, "eta": eta, "x": x, "y": y, "z": z, "a": b["a_eta"],
        "calH": -1 / eta, "geom": geom, "Lam": 3 * H**2,
    }


def linearized_dE(s, h):
    return linearized_einstein_covariant(s["geom"], h, s["Lam"])["dE"]
