"""Canonical variables derived from the quadratic action.

From docs/quadratic_action.md the tensor action in the convention used here is

    S2 = (M_pl^2 / 8) int d eta d^3x a^2 [ gamma_ij' gamma_ij' - d_k gamma_ij d_k gamma_ij ].

With gamma_ij = sum_lam e^lam_ij gamma_lam and e^lam_ij e^lam'_ij = 2 delta,

    S2 = sum_lam (1/2) int d eta d^3x z_T^2 [ gamma_lam'^2 - (d gamma_lam)^2 ],   z_T^2 = M_pl^2 a^2 / 2.

Momentum  pi_lam = z_T^2 gamma_lam'.   Canonical variable v = z_T gamma_lam gives
(after one integration by parts)

    S2 = sum_lam (1/2) int d eta d^3x [ v'^2 - (d v)^2 + (z_T''/z_T) v^2 ],   z_T''/z_T = a''/a = 2/eta^2.

The prefactor is *derived* here from the symbolic Lagrangian produced by the
second-order expansion of the Einstein-Hilbert action, not quoted.
"""

from __future__ import annotations

import functools

import sympy as sp

from ..symbolic import quadratic_action as qa


@functools.lru_cache(maxsize=None)
def derive_z_T_squared():
    """Extract z_T^2 from the verified tensor Lagrangian of ``symbolic.quadratic_action``.

    Returns (z2, details): z2 is the coefficient such that
    (1/(16 pi G)) * L2_target = (1/2) z2 (f'^2 - (d f)^2) for one polarization f.
    """
    s = qa.setup()
    eta, a, fp, Z = s["eta"], s["a"], s["fp"], s["z"]
    M = sp.Symbol("M_pl", positive=True)
    L_target = qa.target_tensor_lagrangian(s)  # coefficient of 1/(16 pi G)
    pref = M**2 / 2                             # 1/(16 pi G) = M_pl^2/2
    L = pref * L_target
    # keep only plus polarization (cross = 0): L = (1/2) z2 (fp'^2 - fp_z^2)
    fc = s["fc"]
    L_plus = L.subs(fc, 0).doit()
    kin = sp.diff(L_plus, sp.diff(fp, eta))  # = z2 * fp'
    z2 = sp.simplify(kin / sp.diff(fp, eta))
    return z2, {"setup": s, "M": M, "L_plus": sp.simplify(L_plus)}


def effective_potential():
    """z''/z for z = a M/sqrt(2) in de Sitter: 2/eta^2."""
    eta = sp.Symbol("eta", negative=True)
    H, M = sp.symbols("H M_pl", positive=True)
    z = (-1 / (H * eta)) * M / sp.sqrt(2)
    return sp.simplify(sp.diff(z, eta, 2) / z)


@functools.lru_cache(maxsize=None)
def verify_canonical_action():
    """Check that L_gamma(gamma = v/z) and L_v differ by a total derivative.

    Uses the Euler-Lagrange operator: if the difference has identically vanishing
    E-L expression with respect to v, it is a total derivative.
    """
    eta = sp.Symbol("eta", negative=True)
    zc = sp.Symbol("z_coord", real=True)
    H, M = sp.symbols("H M_pl", positive=True)
    v = sp.Function("v")(eta, zc)
    z = (-1 / (H * eta)) * M / sp.sqrt(2)
    gam = v / z
    L_gamma = sp.Rational(1, 2) * z**2 * (sp.diff(gam, eta) ** 2 - sp.diff(gam, zc) ** 2)
    zpp_over_z = sp.simplify(sp.diff(z, eta, 2) / z)
    L_v = sp.Rational(1, 2) * (sp.diff(v, eta) ** 2 - sp.diff(v, zc) ** 2 + zpp_over_z * v**2)
    diff = sp.expand(L_gamma - L_v)
    el = qa.euler_expression(sp.simplify(diff), v)
    return sp.simplify(el), zpp_over_z


def mode_equation_from_action():
    """E-L equation of L_v: v'' - v_zz - (z''/z) v = 0  (Fourier: v'' + (k^2 - 2/eta^2) v = 0)."""
    eta = sp.Symbol("eta", negative=True)
    zc = sp.Symbol("z_coord", real=True)
    v = sp.Function("v")(eta, zc)
    zpp = effective_potential()
    L_v = sp.Rational(1, 2) * (sp.diff(v, eta) ** 2 - sp.diff(v, zc) ** 2 + zpp * v**2)
    return sp.simplify(qa.euler_expression(L_v, v)), eta, v, zc
