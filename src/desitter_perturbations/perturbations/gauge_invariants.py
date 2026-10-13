"""Gauge-invariant combinations (flat de Sitter, conformal time).

With sigma = E' - B the Bardeen potentials and the vector invariant are

    Phi_B = A - (sigma' + calH sigma),       Psi_B = psi + calH sigma,
    V_i   = S_i - F_i',                      gamma_ij  (invariant by itself).

Their invariance follows from the transformation laws in
``symbolic.gauge_transformations`` (verified symbolically there) and is verified
again numerically in ``numeric_gauge_experiment``.

Caveat: these are invariant under linearized diffeomorphisms of the *background*
and are therefore well-defined perturbative quantities; this does not make them
local observables in a gravitational theory (see docs/limitations.md).
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def bardeen_symbolic(A, B, psi, E, eta):
    cH = -1 / eta
    sig = sp.diff(E, eta) - B
    Phi = A - (sp.diff(sig, eta) + cH * sig)
    Psi = psi + cH * sig
    return Phi, Psi


def vector_invariant_symbolic(S, F, eta):
    return S - sp.diff(F, eta)


# ----------------------------------------------------- numeric gauge experiment

def gauge_shift_fields(fields, T, L, dT, dL, *, eta):
    """Apply the finite-element law  (A,B,psi,E) -> shifted by gauge functions.

    ``fields`` is a dict of arrays on the eta grid: A, B, psi, E and their
    first eta-derivative E_p (needed for sigma').  ``T, L`` are the gauge functions
    on the grid and ``dT, dL`` their eta-derivatives.  Returns new fields.
    """
    cH = -1.0 / eta
    out = dict(fields)
    out["A"] = fields["A"] - (cH * T + dT)
    out["B"] = fields["B"] + T - dL
    out["psi"] = fields["psi"] + cH * T
    out["E"] = fields["E"] - L
    return out


def bardeen_numeric(fields, eta, dE_deta, d2E_deta2, dB_deta):
    """Numerical Phi_B, Psi_B from arrays (A, B, psi, E and needed derivatives)."""
    cH = -1.0 / eta
    sigma = dE_deta - fields["B"]
    dsigma = d2E_deta2 - dB_deta
    Phi = fields["A"] - (dsigma + cH * sigma)
    Psi = fields["psi"] + cH * sigma
    return Phi, Psi


def numeric_gauge_experiment(n=801, eta_range=(-50.0, -0.05), k=1.0, seed=3):
    """Build one physical scalar configuration, move it through several gauges.

    A smooth reference configuration (A,B,psi,E as analytic functions of eta, times
    the common Fourier factor cos(kz)) is shifted by analytic gauge functions T, L
    (derivatives known exactly).  Returns the coordinate-dependent fields in each
    gauge and the Bardeen potentials, which must coincide across gauges.

    The configuration is *arbitrary* (not a vacuum solution): the gauge-invariance of
    Phi_B, Psi_B is a kinematic statement independent of the field equations.
    """
    eta = np.linspace(eta_range[0], eta_range[1], n)
    # analytic reference fields and derivatives
    A = np.sin(0.3 * eta) / (1 + eta**2 / 25)
    B = np.cos(0.2 * eta) * 0.5
    psi = np.exp(0.02 * eta) * 0.3
    E = 0.1 * np.sin(0.15 * eta)
    dE = 0.1 * 0.15 * np.cos(0.15 * eta)
    d2E = -0.1 * 0.15**2 * np.sin(0.15 * eta)
    dB = -0.5 * 0.2 * np.sin(0.2 * eta)
    # explicit analytic derivatives for each spec (T, L, dT, dL, d2T, d2L)
    def spec_funcs(name):
        if name == "original":
            z0 = 0 * eta
            return z0, z0, z0, z0, z0, z0
        if name == "shifted_1":
            T = 0.4 * np.sin(0.1 * eta); dT = 0.04 * np.cos(0.1 * eta); d2T = -0.004 * np.sin(0.1 * eta)
            L = 0.3 * np.cos(0.1 * eta); dL = -0.03 * np.sin(0.1 * eta); d2L = -0.003 * np.cos(0.1 * eta)
            return T, L, dT, dL, d2T, d2L
        if name == "shifted_2":
            T = 0.2 * np.exp(0.05 * eta); dT = 0.2 * 0.05 * np.exp(0.05 * eta); d2T = 0.2 * 0.05**2 * np.exp(0.05 * eta)
            L = 0.1 * np.sin(0.07 * eta); dL = 0.1 * 0.07 * np.cos(0.07 * eta); d2L = -0.1 * 0.07**2 * np.sin(0.07 * eta)
            return T, L, dT, dL, d2T, d2L
        raise KeyError(name)

    results = {}
    for name in ("original", "shifted_1", "shifted_2"):
        T, L, dT, dL, d2T, d2L = spec_funcs(name)
        cH = -1.0 / eta
        # fields and the derivatives needed for Bardeen
        A_g = A - (cH * T + dT)
        B_g = B + T - dL
        psi_g = psi + cH * T
        E_g = E - L
        dE_g = dE - dL
        d2E_g = d2E - d2L
        dB_g = dB + dT - d2L
        sigma = dE_g - B_g
        dsigma = d2E_g - dB_g
        Phi = A_g - (dsigma + cH * sigma)
        Psi = psi_g + cH * sigma
        results[name] = dict(A=A_g, B=B_g, psi=psi_g, E=E_g, Phi=Phi, Psi=Psi, T=T, L=L)
    return eta, results
