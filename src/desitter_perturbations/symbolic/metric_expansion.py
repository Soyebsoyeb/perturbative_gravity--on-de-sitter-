"""Perturbative expansion of the inverse metric and of the determinant.

Results (indices of h raised with the background metric, h = gbar^{mn} h_{mn}):

    g^{mn}   = gbar^{mn} - eps h^{mn} + eps^2 h^m_r h^{rn} + O(eps^3)
    sqrt(-g) = sqrt(-gbar) [1 + eps h/2 + eps^2 (h^2/8 - h_{mn} h^{mn}/4)] + O(eps^3)

Both are verified here with *generic* symbolic perturbations (ten independent
symbols for h) on a generic diagonal background, so the checks are not tied to
a special metric.  A deliberately wrong variant of each formula is provided so
the test-suite can demonstrate that the checks detect index-placement, factor-
of-two and sign errors (mutation testing).
"""

from __future__ import annotations

import itertools
import sympy as sp

EPS = sp.Symbol("epsilon")


def generic_perturbation(n: int = 4, name: str = "h"):
    """Symmetric matrix of independent symbols h_{ij}."""
    H = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            s = sp.Symbol(f"{name}{i}{j}", real=True)
            H[i, j] = s
            H[j, i] = s
    return H


def generic_diagonal_background(n: int = 4, signature=(-1, 1, 1, 1)):
    """Diagonal background with independent positive functions (symbols) a_i."""
    syms = [sp.Symbol(f"b{i}", positive=True) for i in range(n)]
    return sp.diag(*[signature[i] * syms[i] for i in range(n)])


def truncate(expr, order: int, eps=EPS):
    """Drop terms of order > ``order`` in eps."""
    expr = sp.expand(expr)
    return sum(expr.coeff(eps, k) * eps**k for k in range(order + 1))


def inverse_metric_series(gbar: sp.Matrix, h: sp.Matrix, order: int = 2, eps=EPS):
    """Truncated Neumann series g^{-1} = sum_k (-eps)^k (gbar^{-1} h)^k gbar^{-1}."""
    ginv = gbar.inv()
    M = ginv * h  # mixed tensor h^m_n
    term = ginv
    total = ginv
    for k in range(1, order + 1):
        term = -eps * (M * term)
        total = total + term
    return total.applyfunc(sp.expand)


def inverse_metric_series_wrong_sign(gbar, h, order=2, eps=EPS):
    """Mutation: second-order term with the wrong sign."""
    ginv = gbar.inv()
    M = ginv * h
    return (ginv - eps * M * ginv - eps**2 * M * M * ginv).applyfunc(sp.expand)


def inverse_metric_series_wrong_index(gbar, h, order=2, eps=EPS):
    """Mutation: contracts the wrong index pair at second order (h_{mr} h^r_n
    instead of h^m_r h^{rn}); differs for non-diagonal backgrounds, and here
    inserts a stray transpose that breaks the identity for generic h."""
    ginv = gbar.inv()
    M = ginv * h
    return (ginv - eps * M * ginv + eps**2 * (M.T * M.T) * ginv).applyfunc(sp.expand)


def sqrt_det_series_coefficients(gbar: sp.Matrix, h: sp.Matrix):
    """Return (c1, c2) with sqrt(-g)/sqrt(-gbar) = 1 + eps c1 + eps^2 c2 + O(eps^3).

    Computed from the exact determinant (independent of the trace formula):
    det(1 + eps M) = 1 + eps tr M + eps^2 (tr(M)^2 - tr(M^2))/2 + ...
    expanded here directly from the characteristic polynomial of M.
    """
    M = gbar.inv() * h
    # Exact: ratio of determinants as a polynomial in eps.
    ratio = sp.expand((sp.eye(gbar.shape[0]) + EPS * M).det())
    c1_det = ratio.coeff(EPS, 1)
    c2_det = ratio.coeff(EPS, 2)
    # sqrt(1 + a eps + b eps^2) = 1 + a eps/2 + (b/2 - a^2/8) eps^2
    c1 = sp.simplify(c1_det / 2)
    c2 = sp.simplify(c2_det / 2 - c1_det**2 / 8)
    return c1, c2


def sqrt_det_series_formula(gbar: sp.Matrix, h: sp.Matrix):
    """Closed form quoted in the text: (h/2, h^2/8 - h_{mn}h^{mn}/4)."""
    M = gbar.inv() * h
    tr = M.trace()
    tr2 = (M * M).trace()  # = h^m_n h^n_m = h_{mn} h^{mn}
    return tr / 2, tr**2 / 8 - tr2 / 4


def sqrt_det_series_formula_wrong_factor(gbar, h):
    """Mutation: factor-of-two error in the h_{mn}h^{mn} term."""
    M = gbar.inv() * h
    tr = M.trace()
    tr2 = (M * M).trace()
    return tr / 2, tr**2 / 8 - tr2 / 2


def check_inverse_metric(gbar, h, order=2, builder=inverse_metric_series):
    """Return the residual matrix (g * g^{-1}_trunc - 1) truncated to ``order``.

    A correct expansion has a residual that is O(eps^{order+1}), so after
    truncation to ``order`` it must vanish identically.
    """
    ginv_t = builder(gbar, h, order)
    prod = (gbar + EPS * h) * ginv_t
    res = (prod - sp.eye(gbar.shape[0])).applyfunc(lambda e: truncate(e, order))
    return res.applyfunc(sp.simplify)


def check_sqrt_det(gbar, h, formula=sqrt_det_series_formula):
    c1, c2 = sqrt_det_series_coefficients(gbar, h)
    f1, f2 = formula(gbar, h)
    return sp.simplify(c1 - f1), sp.simplify(c2 - f2)


def two_dimensional_example():
    """Explicit 2x2 example with a non-trivial background, as an independent test.

    gbar = diag(-1, a^2), h arbitrary: confirms the 4D formulas in a setting
    that can be done by hand.
    """
    a = sp.Symbol("a", positive=True)
    gbar = sp.diag(-1, a**2)
    h = generic_perturbation(2)
    return gbar, h
