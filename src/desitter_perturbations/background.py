"""Exact de Sitter background in flat slicing: cosmic-time and conformal-time forms.

Relation between the two time variables
---------------------------------------
With a(t) = exp(H t) and d eta = dt / a:

    eta = int dt exp(-H t) = -exp(-H t) / H   (integration constant chosen so eta -> 0^- as t -> inf)

so that a(eta) = -1/(H eta), eta < 0.  The conformal Hubble parameter is
calH = a'/a = -1/eta = a H.  The comoving Hubble radius is 1/(a H) = -eta.
"""

from __future__ import annotations

import numpy as np
import sympy as sp

from .conventions import Units, DEFAULT_UNITS

# ---------------------------------------------------------------- numerical

def scale_factor_t(t, H: float = 1.0):
    return np.exp(H * np.asarray(t, dtype=float))


def conformal_time_from_t(t, H: float = 1.0):
    """eta(t) = -exp(-H t)/H."""
    return -np.exp(-H * np.asarray(t, dtype=float)) / H


def t_from_conformal_time(eta, H: float = 1.0):
    eta = np.asarray(eta, dtype=float)
    if np.any(eta >= 0):
        raise ValueError("conformal time must be negative in flat de Sitter slicing")
    return -np.log(-H * eta) / H


def scale_factor_eta(eta, H: float = 1.0):
    eta = np.asarray(eta, dtype=float)
    return -1.0 / (H * eta)


def conformal_hubble(eta):
    """calH = a'/a = -1/eta."""
    return -1.0 / np.asarray(eta, dtype=float)


def x_of_eta(k, eta):
    """Horizon variable x = k/(aH) = -k eta."""
    return -np.asarray(k, dtype=float) * np.asarray(eta, dtype=float)


def eta_of_x(k, x):
    return -np.asarray(x, dtype=float) / np.asarray(k, dtype=float)


def physical_wavelength(k, eta, H: float = 1.0):
    """2 pi a / k."""
    return 2.0 * np.pi * scale_factor_eta(eta, H) / k


def hubble_radius_comoving(eta):
    """1/(aH) = -eta."""
    return -np.asarray(eta, dtype=float)


# ----------------------------------------------------------------- symbolic

t, eta, H = sp.symbols("t eta H", real=True)
x, y, z = sp.symbols("x y z", real=True)
H_pos = sp.Symbol("H", positive=True)


def symbolic_background(H_sym=H_pos):
    """Return dictionary of symbolic objects describing both coordinate systems."""
    t_s = sp.Symbol("t", real=True)
    eta_s = sp.Symbol("eta", negative=True)
    X = sp.symbols("x y z", real=True)
    a_t = sp.exp(H_sym * t_s)
    a_eta = -1 / (H_sym * eta_s)
    g_cosmic = sp.diag(-1, a_t**2, a_t**2, a_t**2)
    g_conf = sp.diag(-a_eta**2, a_eta**2, a_eta**2, a_eta**2)
    return {
        "H": H_sym,
        "t": t_s,
        "eta": eta_s,
        "spatial": X,
        "a_t": a_t,
        "a_eta": a_eta,
        "coords_cosmic": (t_s,) + tuple(X),
        "coords_conformal": (eta_s,) + tuple(X),
        "g_cosmic": g_cosmic,
        "g_conformal": g_conf,
        "calH": sp.simplify(sp.diff(a_eta, eta_s) / a_eta),
    }


def verify_time_relation(H_sym=H_pos):
    """Symbolic check of d eta/dt = 1/a and a(eta(t)) = exp(H t)."""
    t_s = sp.Symbol("t", real=True)
    eta_of_t = -sp.exp(-H_sym * t_s) / H_sym
    a_t = sp.exp(H_sym * t_s)
    r1 = sp.simplify(sp.diff(eta_of_t, t_s) - 1 / a_t)
    r2 = sp.simplify((-1 / (H_sym * eta_of_t)) - a_t)
    return r1, r2
