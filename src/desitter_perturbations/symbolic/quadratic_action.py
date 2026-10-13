"""Second-order Einstein-Hilbert action for transverse-traceless tensor modes.

Setup
-----
Background: flat de Sitter in conformal time, a = -1/(H eta).
Perturbation (exact in eps, truncated afterwards):

    ds^2 = a^2 [ -d eta^2 + (delta_ij + eps gamma_ij) dx^i dx^j ],
    gamma_xx = -gamma_yy = f_plus(eta, z),  gamma_xy = f_cross(eta, z),

i.e. a wave travelling along z with both polarizations, using e^+_{ij} = diag(1,-1,0),
e^x_{ij} (off-diagonal), normalised so that e^lam_ij e^lam_ij = 2.  Then
gamma_ij gamma_ij = 2 (f_plus^2 + f_cross^2).

The Lagrangian density (coefficient of 1/(16 pi G)) is  L = sqrt(-g) (R - 2 Lambda).

Route 1 (direct expansion): expand L exactly to eps^2.
Route 2 (linearized equations): L2 = -(1/2) sqrt(-gbar) h^{mn} dE_{mn}[h], valid on a
    background solution because  delta S = -(1/16 pi G) int sqrt(-g) E^{mn} delta g_{mn}
    and Ebar = 0.

Target (derived below and *verified* here, not assumed):

    L2 = (a^2/4) (gamma_ij' gamma_ij' - d_k gamma_ij d_k gamma_ij) + total derivative,

so S2 = (1/(16 pi G)) int L2 = (M_pl^2/8) int d eta d^3x a^2 [gamma_ij'^2 - (d gamma_ij)^2].
"""

from __future__ import annotations

import functools

import sympy as sp

from .curvature import Geometry
from .einstein_tensor import linearized_einstein_covariant
from .metric_expansion import EPS


def setup():
    eta = sp.Symbol("eta", negative=True)
    H = sp.Symbol("H", positive=True)
    X, Y, Z = sp.symbols("x y z", real=True)
    coords = (eta, X, Y, Z)
    a = -1 / (H * eta)
    fp = sp.Function("f_p")(eta, Z)
    fc = sp.Function("f_c")(eta, Z)
    gam = sp.Matrix([[0, 0, 0, 0], [0, fp, fc, 0], [0, fc, -fp, 0], [0, 0, 0, 0]])
    gbar = sp.diag(-a**2, a**2, a**2, a**2)
    h = a**2 * gam
    return dict(eta=eta, H=H, coords=coords, a=a, fp=fp, fc=fc, gam=gam, gbar=gbar, h=h, z=Z)


def _inverse_exact(g):
    """Exact inverse for the block structure used here (sympy handles 4x4 fine)."""
    return sp.simplify(g.inv())


@functools.lru_cache(maxsize=None)
def route1_direct_expansion(s=None):
    """Return L2 = coefficient of eps^2 in sqrt(-g)(R - 2 Lambda), plus L0, L1."""
    s = s or setup()
    eta, H, coords = s["eta"], s["H"], s["coords"]
    g = s["gbar"] + EPS * s["h"]
    ginv = _inverse_exact(g)
    geom = Geometry(coords, g, ginv=ginv, simp=lambda e: sp.cancel(sp.together(e)))
    R = geom.ricci_scalar
    detg = sp.simplify(g.det())
    sqrtg = sp.sqrt(sp.factor(-detg))
    Lam = 3 * H**2
    L = sqrtg * (R - 2 * Lam)
    # Taylor coefficients at eps = 0
    L0 = sp.simplify(L.subs(EPS, 0))
    L1 = sp.simplify(sp.diff(L, EPS).subs(EPS, 0))
    L2 = sp.simplify(sp.diff(L, EPS, 2).subs(EPS, 0) / 2)
    return {"L0": L0, "L1": L1, "L2": L2, "setup": s}


def target_tensor_lagrangian(s=None):
    """(a^2/4)(gamma_ij' gamma_ij' - d_k gamma_ij d_k gamma_ij) for the ansatz."""
    s = s or setup()
    eta, a, fp, fc, Z = s["eta"], s["a"], s["fp"], s["fc"], s["z"]
    gsq_t = 2 * (sp.diff(fp, eta)**2 + sp.diff(fc, eta)**2)
    gsq_z = 2 * (sp.diff(fp, Z)**2 + sp.diff(fc, Z)**2)
    return a**2 / 4 * (gsq_t - gsq_z)


@functools.lru_cache(maxsize=None)
def route2_from_linearized_equations(s=None):
    """L2 = -(1/2) sqrt(-gbar) h^{mn} dE_{mn}[h] using the covariant linearization."""
    s = s or setup()
    eta, H, coords = s["eta"], s["H"], s["coords"]
    geom = Geometry(coords, s["gbar"])
    out = linearized_einstein_covariant(geom, s["h"], 3 * H**2)
    dE = out["dE"]
    ginv = geom.ginv
    h_up = ginv * s["h"] * ginv
    sqrt_gbar = s["a"]**4
    L2 = -sp.Rational(1, 2) * sqrt_gbar * sum(h_up[m, n] * dE[m, n] for m in range(4) for n in range(4))
    return {"L2": sp.simplify(L2), "dE": dE, "setup": s}


def euler_expression(L, f):
    """Euler-Lagrange expression  dL/df - D_a dL/d(d_a f) + D_a D_b dL/d(d_a d_b f) - ...

    Implemented directly (rather than via sympy.calculus.euler, which silently
    returns an empty list for some one-variable cases) and valid for any order.
    """
    res = sp.diff(L, f)
    derivs = [d for d in L.atoms(sp.Derivative) if d.expr == f]
    for d in derivs:
        vars_ = []
        for v, n in d.variable_count:
            vars_ += [v] * n
        term = sp.diff(L, d)
        for v in vars_:
            term = sp.diff(term, v)
        res += (-1) ** len(vars_) * term
    return res


def euler_operator(L, funcs, coords_eta_z=None):
    """Euler-Lagrange expressions of L with respect to each function in ``funcs``."""
    return [sp.simplify(euler_expression(L, f)) for f in funcs]


def total_derivative_check(Ldiff, s=None):
    """A Lagrangian is a total derivative iff its Euler-Lagrange expressions vanish
    identically.  Returns the list of E-L expressions of ``Ldiff``."""
    s = s or setup()
    return euler_operator(Ldiff, [s["fp"], s["fc"]], [s["eta"], s["z"]])


def tensor_equation_from_action(s=None):
    """E-L equations of the target Lagrangian; returns (eq_plus, eq_cross)."""
    s = s or setup()
    return euler_operator(target_tensor_lagrangian(s), [s["fp"], s["fc"]], [s["eta"], s["z"]])


@functools.lru_cache(maxsize=None)
def first_order_variation_check():
    """Bulk first-order variation vanishes on the background solution.

    For a diagonal perturbation h_{mn} = a^2 diag(-2 A, c1, c2, c3)(eta) with arbitrary
    functions, the O(eps) part of sqrt(-g)(R - 2 Lambda) must be a total derivative
    (its Euler-Lagrange expressions vanish identically), because
    L1 = -sqrt(-gbar) Ebar^{mn} h_{mn} + total derivative and Ebar = 0.
    Also returns the same coefficient for a deliberately *wrong* Lambda (Lambda = 2 H^2)
    to show the check has teeth: the E-L expressions are then non-zero.
    """
    eta = sp.Symbol("eta", negative=True)
    H = sp.Symbol("H", positive=True)
    X, Y, Z = sp.symbols("x y z", real=True)
    coords = (eta, X, Y, Z)
    a = -1 / (H * eta)
    A = sp.Function("A")(eta)
    c = [sp.Function(f"c{i}")(eta) for i in (1, 2, 3)]
    gbar = sp.diag(-a**2, a**2, a**2, a**2)
    h = a**2 * sp.diag(-2 * A, *c)
    g = gbar + EPS * h
    ginv = sp.simplify(g.inv())
    geom = Geometry(coords, g, ginv=ginv, simp=lambda e: sp.cancel(sp.together(e)))
    R = geom.ricci_scalar
    sqrtg = sp.sqrt(sp.factor(-sp.simplify(g.det())))
    out = {}
    for label, Lam in (("correct", 3 * H**2), ("wrong", 2 * H**2)):
        L = sqrtg * (R - 2 * Lam)
        L1 = sp.simplify(sp.diff(L, EPS).subs(EPS, 0))
        out[label] = {"L1": L1, "EL": euler_operator(L1, [A] + c)}
    return out
