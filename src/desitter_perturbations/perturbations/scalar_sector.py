"""Scalar sector in a general gauge and its reduction to Bardeen potentials.

Metric (fields depend on (eta, z); all gauge functions along z are enough to
expose the structure because the general case follows by rotational symmetry):

    h_00 = -2 a^2 A,   h_0z = a^2 d_z B,   h_xx = h_yy = -2 a^2 psi,
    h_zz = a^2 (-2 psi + 2 d_z^2 E).

Verified results (exact symbolic identities, see tests/test_svt_decomposition.py):

    dE_00 = 2 nabla^2 Psi_B - 6 calH (Psi_B' + calH Phi_B)
    dE_0z = 2 d_z (Psi_B' + calH Phi_B)
    dE_xx - dE_zz = d_z^2 (Phi_B - Psi_B)

so every component of the linearized field equation depends on the metric
perturbation only through the gauge invariants Phi_B and Psi_B.
"""

from __future__ import annotations

import functools

import sympy as sp

from ._common import setup, linearized_dE
from .gauge_invariants import bardeen_symbolic


@functools.lru_cache(maxsize=None)
def general_gauge_dE():
    s = setup()
    eta, z, a = s["eta"], s["z"], s["a"]
    A = sp.Function("A")(eta, z)
    B = sp.Function("B")(eta, z)
    psi = sp.Function("psi")(eta, z)
    E = sp.Function("E")(eta, z)
    h = sp.zeros(4, 4)
    h[0, 0] = -2 * a**2 * A
    h[0, 3] = h[3, 0] = a**2 * sp.diff(B, z)
    h[1, 1] = h[2, 2] = -2 * a**2 * psi
    h[3, 3] = a**2 * (-2 * psi + 2 * sp.diff(E, z, 2))
    dE = linearized_dE(s, h)
    Phi, Psi = bardeen_symbolic(A, B, psi, E, eta)
    cH = s["calH"]
    pred = {
        (0, 0): 2 * sp.diff(Psi, z, 2) - 6 * cH * (sp.diff(Psi, eta) + cH * Phi),
        (0, 3): 2 * sp.diff(sp.diff(Psi, eta) + cH * Phi, z),
    }
    traceless_pred = sp.diff(Phi - Psi, z, 2)
    residuals = {
        "00": sp.simplify(dE[0, 0] - pred[(0, 0)]),
        "0z": sp.simplify(dE[0, 3] - pred[(0, 3)]),
        "xx-zz": sp.simplify(dE[1, 1] - dE[3, 3] - traceless_pred),
    }
    return {"dE": dE, "Phi": Phi, "Psi": Psi, "residuals": residuals, "setup": s,
            "fields": (A, B, psi, E)}


@functools.lru_cache(maxsize=None)
def newtonian_gauge_dE():
    """dE components in Newtonian gauge (B = E = 0): A -> Phi, psi -> Psi."""
    s = setup()
    eta, z, a = s["eta"], s["z"], s["a"]
    Phi = sp.Function("Phi")(eta, z)
    Psi = sp.Function("Psi")(eta, z)
    h = a**2 * sp.diag(-2 * Phi, -2 * Psi, -2 * Psi, -2 * Psi)
    dE = linearized_dE(s, h)
    return {"dE": dE, "Phi": Phi, "Psi": Psi, "setup": s}
