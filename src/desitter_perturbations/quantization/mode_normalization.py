"""Mode normalisation from the canonical commutation relation (symbolic derivation)."""

from __future__ import annotations

import sympy as sp


def normalization_from_wronskian():
    """Fix |C|^2 in gamma = C (1 + i k eta) exp(-i k eta) from z_T^2 (gamma gamma*' - c.c.) = i.

    Returns (C_abs_squared, residual).  The result H^2/(M_pl^2 k^3) is the origin of the
    amplitude in ``dynamics.analytic_solutions``.
    """
    k, H, M = sp.symbols("k H M_pl", positive=True)
    eta = sp.Symbol("eta", negative=True)
    C = sp.Symbol("C", positive=True)  # modulus; phase is unobservable
    g = C * (1 + sp.I * k * eta) * sp.exp(-sp.I * k * eta)
    gs = C * (1 - sp.I * k * eta) * sp.exp(sp.I * k * eta)
    z2 = M**2 / (2 * H**2 * eta**2)
    W = sp.simplify(z2 * (g * sp.diff(gs, eta) - gs * sp.diff(g, eta)))
    sol = sp.solve(sp.Eq(W, sp.I), C**2)
    C2 = sp.simplify(sol[0]) if sol else None
    return C2, W


def verify_mode_solves_equation():
    """Residual of gamma'' - (2/eta) gamma' + k^2 gamma for gamma = (1 + i k eta) e^{-i k eta}."""
    k = sp.Symbol("k", positive=True)
    eta = sp.Symbol("eta", negative=True)
    g = (1 + sp.I * k * eta) * sp.exp(-sp.I * k * eta)
    return sp.simplify(sp.diff(g, eta, 2) - 2 / eta * sp.diff(g, eta) + k**2 * g)


def positive_frequency_limit():
    """lim_{k eta -> -inf} sqrt(2k) v_BD e^{i k eta} = 1 (BD selects the e^{-i k eta} branch)."""
    k = sp.Symbol("k", positive=True)
    eta = sp.Symbol("eta", negative=True)
    v = sp.exp(-sp.I * k * eta) / sp.sqrt(2 * k) * (1 - sp.I / (k * eta))
    expr = sp.simplify(sp.sqrt(2 * k) * v * sp.exp(sp.I * k * eta))
    return sp.limit(expr, eta, -sp.oo)
