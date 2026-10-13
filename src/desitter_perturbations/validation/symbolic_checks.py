"""Exact (symbolic) identity checks.  Each function returns a list of Records.

Residuals of symbolic checks are 0.0 (identically zero after simplification) or 1.0.
Negative controls demonstrate that the checks are sensitive to the errors they
claim to detect (wrong index placement, factor of two, curvature sign, wrong Lambda).
"""

from __future__ import annotations

import functools
import sympy as sp

from .records import record, zero_residual
from .. import background as bg
from ..symbolic import (metric_expansion as me, curvature as cu, einstein_tensor as et,
                        bianchi_identities as bi, gauge_transformations as gt,
                        quadratic_action as qa, tensor_decomposition as td)
from ..perturbations import scalar_sector, vector_sector, tensor_sector, constraints
from ..quantization import canonical_variables as cv, mode_normalization as mn


def _bg():
    return bg.symbolic_background()


@functools.lru_cache(maxsize=None)
def _geometries():
    b = _bg()
    gc = cu.Geometry(b["coords_cosmic"], b["g_cosmic"])
    gf = cu.Geometry(b["coords_conformal"], b["g_conformal"])
    return b, gc, gf


@functools.lru_cache(maxsize=None)
def _generic_linearization():
    """Covariant and brute-force dE for a fixed generic polynomial perturbation (cached)."""
    b, _, gf = _geometries()
    eta = b["eta"]
    x, y, z = b["spatial"]
    H = b["H"]
    h = sp.Matrix([[eta * x, y, 0, z * eta], [y, x**2, eta * z, 0],
                   [0, eta * z, y * z, x], [z * eta, 0, x, eta**2]])
    Lam = 3 * H**2
    cov = et.linearized_einstein_covariant(gf, h, Lam)
    e0, e1 = et.linearized_einstein_bruteforce(b["coords_conformal"], b["g_conformal"], h, Lam)
    return h, cov, e0, e1


def background_checks():
    recs = []
    b, gc, gf = _geometries()
    H = b["H"]
    r1, r2 = bg.verify_time_relation()
    recs.append(record("bg.time_relation", "background",
                       "d eta/dt = 1/a and a(eta(t)) = exp(H t) with eta = -exp(-H t)/H",
                       zero_residual([r1, r2]), 0, method="sympy simplify"))
    for name, G in (("cosmic", gc), ("conformal", gf)):
        R = G.ricci_scalar
        recs.append(record(f"bg.ricci_scalar.{name}", "background", f"Ricci scalar = 12 H^2 ({name} slicing)",
                           zero_residual(R - 12 * H**2), 0, method="direct Christoffel/Riemann contraction"))
        recs.append(record(f"bg.ricci_tensor.{name}", "background", f"R_mn = 3 H^2 g_mn ({name})",
                           zero_residual(G.ricci - 3 * H**2 * G.g), 0))
        Gt = et.einstein_tensor(G)
        recs.append(record(f"bg.einstein.{name}", "background", f"G_mn = -3 H^2 g_mn ({name})",
                           zero_residual(Gt + 3 * H**2 * G.g), 0))
        recs.append(record(f"bg.field_eq.{name}", "background", f"G_mn + Lambda g_mn = 0 with Lambda = 3H^2 ({name})",
                           zero_residual(et.field_equation(G, 3 * H**2)), 0))
        res = G.constant_curvature_residual(H)
        recs.append(record(f"bg.constant_curvature.{name}", "background",
                           f"R_mnrs = H^2 (g_mr g_ns - g_ms g_nr), all 256 components ({name})",
                           zero_residual(res), 0, method="component comparison"))
        recs.append(record(f"bg.kretschmann.{name}", "background", f"R_abcd R^abcd = 24 H^4 ({name})",
                           zero_residual(G.kretschmann() - 24 * H**4), 0))
        recs.append(record(f"bg.ricci_squared.{name}", "background", f"R_ab R^ab = 36 H^4 ({name})",
                           zero_residual(G.ricci_squared() - 36 * H**4), 0))
    # sign negative control
    Gneg = cu.Geometry(b["coords_conformal"], b["g_conformal"], sign=-1)
    recs.append(record("bg.sign_control", "background",
                       "Flipped Riemann sign gives R = -12 H^2 (convention test is sensitive)",
                       zero_residual(Gneg.ricci_scalar - 12 * H**2), 0, expect_fail=True))
    # coordinate-independence: same scalar invariants from two slicings
    recs.append(record("bg.slicing_invariants", "background",
                       "Curvature scalars agree between cosmic-time and conformal-time slicings",
                       zero_residual([gc.ricci_scalar - gf.ricci_scalar, gc.kretschmann() - gf.kretschmann()]), 0))
    # embedding check
    Hs = sp.Symbol("H", positive=True)
    t = sp.Symbol("t", real=True)
    X, Y, Z = sp.symbols("x y z", real=True)
    r2 = X**2 + Y**2 + Z**2
    X0 = sp.sinh(Hs * t) / Hs + Hs / 2 * sp.exp(Hs * t) * r2
    X4 = sp.cosh(Hs * t) / Hs - Hs / 2 * sp.exp(Hs * t) * r2
    emb = [X0, sp.exp(Hs * t) * X, sp.exp(Hs * t) * Y, sp.exp(Hs * t) * Z, X4]
    coords = (t, X, Y, Z)
    sig = [-1, 1, 1, 1, 1]
    ind = sp.zeros(4, 4)
    for m in range(4):
        for n in range(4):
            ind[m, n] = sum(sig[A] * sp.diff(emb[A], coords[m]) * sp.diff(emb[A], coords[n]) for A in range(5))
    target = sp.diag(-1, sp.exp(2 * Hs * t), sp.exp(2 * Hs * t), sp.exp(2 * Hs * t))
    hyper = sum(sig[A] * emb[A] ** 2 for A in range(5))
    recs.append(record("bg.embedding", "background",
                       "Flat slicing embeds in the hyperboloid -X0^2 + X.X + X4^2 = 1/H^2 with induced metric -dt^2 + a^2 dx^2",
                       zero_residual([(ind - target).applyfunc(lambda e: e.rewrite(sp.exp)),
                                          hyper.rewrite(sp.exp) - 1 / Hs**2]), 0,
                       method="pull-back of the 5D Minkowski metric"))
    return recs


def metric_expansion_checks():
    recs = []
    gbar = me.generic_diagonal_background()
    h = me.generic_perturbation()
    res = me.check_inverse_metric(gbar, h, order=2)
    recs.append(record("exp.inverse_metric", "metric_expansion",
                       "g^{mn} = gbar^{mn} - eps h^{mn} + eps^2 h^m_r h^{rn} + O(eps^3) (generic 4D, 10 symbols)",
                       zero_residual(res), 0, assumptions="diagonal background with 4 independent positive symbols; generic symmetric h"))
    res_wrong = me.check_inverse_metric(gbar, h, order=2, builder=me.inverse_metric_series_wrong_sign)
    recs.append(record("exp.inverse_metric.neg_sign", "metric_expansion",
                       "Wrong sign at second order is detected", zero_residual(res_wrong), 0, expect_fail=True))
    # index placement mutation on a background where it matters (non-diagonal)
    gbar_nd = sp.Matrix([[-2, 1, 0, 0], [1, 3, 0, 0], [0, 0, 1, 0], [0, 0, 0, 2]])
    res_idx = me.check_inverse_metric(gbar_nd, h, order=2, builder=me.inverse_metric_series_wrong_index)
    recs.append(record("exp.inverse_metric.neg_index", "metric_expansion",
                       "Transposed index contraction at second order is detected (non-diagonal background)",
                       zero_residual(res_idx), 0, expect_fail=True))
    res_nd = me.check_inverse_metric(gbar_nd, h, order=2)
    recs.append(record("exp.inverse_metric.nondiagonal", "metric_expansion",
                       "Inverse-metric expansion also holds on a non-diagonal background",
                       zero_residual(res_nd), 0))
    d1, d2 = me.check_sqrt_det(gbar, h)
    recs.append(record("exp.sqrt_det", "metric_expansion",
                       "sqrt(-g)/sqrt(-gbar) = 1 + eps h/2 + eps^2 (h^2/8 - h_mn h^mn/4) from the exact determinant",
                       zero_residual([d1, d2]), 0, method="characteristic polynomial of gbar^{-1} h"))
    w1, w2 = me.check_sqrt_det(gbar, h, formula=me.sqrt_det_series_formula_wrong_factor)
    recs.append(record("exp.sqrt_det.neg_factor", "metric_expansion",
                       "Factor-of-two error in the h_mn h^mn term is detected", zero_residual([w1, w2]), 0,
                       expect_fail=True))
    g2, h2 = me.two_dimensional_example()
    res2 = me.check_inverse_metric(g2, h2, order=2)
    e1, e2 = me.check_sqrt_det(g2, h2)
    recs.append(record("exp.two_dim_example", "metric_expansion",
                       "Independent explicit 2D example (gbar = diag(-1, a^2))",
                       zero_residual([res2, e1, e2]), 0))
    return recs


def linearized_einstein_checks():
    recs = []
    b, gc, gf = _geometries()
    h, cov, e0, e1 = _generic_linearization()
    diff = (cov["dE"] - e1).applyfunc(lambda e: sp.cancel(sp.together(e)))
    recs.append(record("lin.background_zero", "linearized_einstein",
                       "Background field equation holds in the truncated brute-force pipeline (O(eps^0) part)",
                       zero_residual(e0), 0))
    recs.append(record("lin.cov_vs_component", "linearized_einstein",
                       "Background-covariant dE_mn equals component-level coefficient of eps in G+Lambda g",
                       zero_residual(diff), 0,
                       assumptions="generic polynomial h (all 10 components non-zero) in conformal coordinates",
                       method="exact rational-function cancellation"))
    # delta Gamma check
    eta = b["eta"]
    H = b["H"]
    gbar = b["g_conformal"]
    ginv_t = me.inverse_metric_series(gbar, h, order=1)
    from ..symbolic.connections import delta_christoffel_direct, delta_christoffel
    dG_direct = delta_christoffel_direct(gbar + me.EPS * h, ginv_t, b["coords_conformal"], me.EPS,
                                         simp=lambda e: sp.cancel(sp.together(e)))
    dG_cov = cov["dGamma"]
    err = [sp.cancel(sp.together(dG_direct[l][m][n] - dG_cov[l][m][n]))
           for l in range(4) for m in range(4) for n in range(4)]
    recs.append(record("lin.delta_gamma", "linearized_einstein",
                       "dGamma^l_mn = 1/2 g^{lr}(nabla_m h_rn + nabla_n h_rm - nabla_r h_mn) vs direct derivative of Gamma",
                       zero_residual(err), 0))
    # Lambda h piece: on shell dE = dG + Lambda h and dR trace relation
    recs.append(record("lin.dE_definition", "linearized_einstein",
                       "dE_mn = dG_mn + Lambda h_mn",
                       zero_residual(cov["dE"] - cov["dG"] - 3 * H**2 * h), 0))
    # Bianchi
    res_cov, naive, corr = bi.bianchi_covariant_residual(gf, cov["dG"], cov["dGamma"])
    recs.append(record("bianchi.covariant_form", "bianchi",
                       "gbar^{mr}[nabla_r dG_mn - dGamma^l_rm Gbar_ln - dGamma^l_rn Gbar_ml] = 0", zero_residual(res_cov), 0))
    recs.append(record("bianchi.naive_divergence_nonzero", "bianchi",
                       "Naive nabla^m dG_mn is NOT zero: connection-variation terms are required",
                       zero_residual(naive), 0, expect_fail=True,
                       notes="the nonzero remainder is proportional to Lambda = -Gbar/g"))
    recs.append(record("bianchi.mixed_form", "bianchi",
                       "Mixed-index linearized Bianchi identity", zero_residual(
                           bi.bianchi_mixed_residual(gf, h, cov["dG"], cov["dGamma"])), 0))
    recs.append(record("bianchi.field_equation_divergence", "bianchi",
                       "nabla^m dE_mn = 0 exactly (no connection-variation terms because Ebar = 0)",
                       zero_residual(bi.bianchi_field_equation_residual(gf, cov["dE"])), 0))
    # cosmic-time background-covariant check on a second coordinate system
    t = b["t"]
    x, y, z = b["spatial"]
    hc = sp.Matrix([[x * t, y, 0, 0], [y, t * z, 0, x], [0, 0, y * t, 0], [0, x, 0, z]])
    covc = et.linearized_einstein_covariant(gc, hc, 3 * H**2)
    e0c, e1c = et.linearized_einstein_bruteforce(b["coords_cosmic"], b["g_cosmic"], hc, 3 * H**2)
    dc = (covc["dE"] - e1c).applyfunc(lambda e: sp.cancel(sp.together(e)))
    recs.append(record("lin.cov_vs_component.cosmic", "linearized_einstein",
                       "Covariant vs component-level linearization in cosmic-time coordinates", zero_residual(dc), 0))
    return recs


def gauge_checks():
    recs = []
    b, gc, gf = _geometries()
    eta = b["eta"]
    x, y, z = b["spatial"]
    H = b["H"]
    xi = [eta**2 * x * y, x * z**2 + eta, y * eta * z, z**3 * x]
    recs.append(record("gauge.lie_equals_covariant", "gauge",
                       "L_xi gbar = nabla_m xi_n + nabla_n xi_m", zero_residual(gt.lie_equals_covariant_residual(xi, gf)), 0))
    hpg = gt.pure_gauge_perturbation(xi, gf)
    dE = et.linearized_einstein_covariant(gf, hpg, 3 * H**2)["dE"]
    recs.append(record("gauge.pure_gauge_solves_linearized", "gauge",
                       "Pure-gauge h = -(nabla xi + nabla xi) satisfies dE_mn = 0", zero_residual(dE), 0))
    # invariance of a non-trivial h: dE[h + h_pg] == dE[h]
    h, cov, _, _ = _generic_linearization()
    dE2 = et.linearized_einstein_covariant(gf, h + hpg, 3 * H**2)["dE"]
    recs.append(record("gauge.field_equation_invariant", "gauge",
                       "dE[h + h_pg] = dE[h] (linearized field equation is gauge invariant)",
                       zero_residual((dE2 - cov["dE"]).applyfunc(lambda e: sp.cancel(sp.together(e)))), 0))
    laws = gt.check_svt_laws()
    recs.append(record("gauge.svt_laws", "gauge",
                       "SVT transformation laws A,psi,B,E,S,F derived from L_xi gbar",
                       zero_residual(laws), 0, method="component matching of -L_xi gbar"))
    # Residual gauge freedom.  Synchronous gauge keeps A = 0 (and B = 0): A -> A - (calH T + T'),
    # so T must satisfy (a T)' = 0, i.e. T = C/a.  Verify that T = C/a leaves A unchanged.
    eta_s = sp.Symbol("eta", negative=True)
    Hs = sp.Symbol("H", positive=True)
    C = sp.Symbol("C")
    a_eta = -1 / (Hs * eta_s)
    T = C / a_eta
    dA = -((-1 / eta_s) * T + sp.diff(T, eta_s))
    recs.append(record("gauge.synchronous_residual", "gauge",
                       "Synchronous gauge keeps A = 0 under T = C/a: residual gauge freedom exists (a T = const)",
                       zero_residual(dA), 0, notes="documented in docs/gauge_symmetry.md; Newtonian and flat gauges have none for k != 0"))
    # Newtonian gauge: B = E = 0 preserved requires L = 0 and T = L' = 0 (k != 0): no residual freedom.
    L = sp.Function("L")(eta_s)
    T2 = sp.Function("T")(eta_s)
    dE_ = -L
    dB_ = T2 - sp.diff(L, eta_s)
    sol = sp.solve([sp.Eq(dE_, 0)], [L], dict=True)[0]
    Tsol = sp.solve(sp.Eq(dB_.subs(L, sol[L]), 0), T2)
    recs.append(record("gauge.newtonian_no_residual", "gauge",
                       "Newtonian gauge (B = E = 0) preserved only by T = L = 0 for k != 0",
                       zero_residual(Tsol[0]) if Tsol else 1.0, 0))
    return recs


def svt_checks():
    recs = []
    res = td.verify_svt_identities()
    for name, flag in td.all_zero(res).items():
        recs.append(record(f"svt.{name}", "svt", f"Fourier-space SVT identity: {name}", 0.0 if flag else 1.0, 0,
                           assumptions="k != 0; generic symmetric h (6 symbols), generic k (3 symbols)"))
    sc = scalar_sector.general_gauge_dE()["residuals"]
    for key, val in sc.items():
        recs.append(record(f"sector.scalar.{key}", "sectors",
                           f"dE component {key} depends on h only through Bardeen potentials", zero_residual(val), 0))
    vc = vector_sector.vector_dE()["residuals"]
    for key, val in vc.items():
        recs.append(record(f"sector.vector.{key}", "sectors",
                           f"Vector-sector dE {key} depends only on V = S - F'", zero_residual(val), 0))
    tr = tensor_sector.tensor_dE()["residuals"]
    recs.append(record("sector.tensor.eom", "sectors",
                       "TT sector: dE_ij = (1/2)(gamma'' + 2 calH gamma' - nabla^2 gamma), other components vanish",
                       zero_residual([tr["xx"], tr["yy"], tr["xy"]] + tr["others"]), 0))
    cte = tensor_sector.covariant_tensor_equation()
    recs.append(record("sector.tensor.covariant", "sectors",
                       "TT ansatz is transverse and traceless, and dE_mn = -(1/2)(box h_mn - 2 H^2 h_mn) (massless spin-2 on de Sitter)",
                       zero_residual([cte["divergence"], cte["trace"], cte["equation"]]), 0,
                       method="two background-covariant derivatives computed componentwise"))
    c = constraints.scalar_constraint_solution()
    recs.append(record("constraints.scalar_traceless", "constraints",
                       "Traceless ij equation gives Phi_B = Psi_B", zero_residual(c["phi_solution"] - c["ps"]), 0))
    k = c["k"]
    recs.append(record("constraints.scalar_vacuum", "constraints",
                       "After 0i and traceless constraints the 00 equation is -2 k^2 Psi_B = 0, so Psi_B = 0 for k != 0",
                       zero_residual(c["e00_after_constraints"] + 2 * k**2 * c["ps"]), 0))
    return recs


def action_checks():
    recs = []
    s = qa.setup()
    r1 = qa.route1_direct_expansion()
    r2 = qa.route2_from_linearized_equations()
    tgt = qa.target_tensor_lagrangian(s)
    recs.append(record("action.background_lagrangian", "action",
                       "Background Lagrangian density sqrt(-g)(R - 2 Lambda) = 6/(H^2 eta^4)",
                       zero_residual(r1["L0"] - 6 / (s["H"] ** 2 * s["eta"] ** 4)), 0))
    recs.append(record("action.first_order_tensor", "action",
                       "First-order term of the TT expansion vanishes identically", zero_residual(r1["L1"]), 0))
    recs.append(record("action.route1_vs_target", "action",
                       "Route 1 (direct expansion): L2 = (a^2/4)(gamma'^2 - (d gamma)^2) + total derivative",
                       zero_residual(qa.total_derivative_check(sp.simplify(r1["L2"] - tgt), s)), 0,
                       method="Euler-Lagrange operator of the difference vanishes <=> total derivative",
                       assumptions="boundary terms discarded: spatial periodicity / compact support, time boundaries ignored"))
    recs.append(record("action.route2_vs_target", "action",
                       "Route 2 (linearized equations): L2 = -(1/2) sqrt(-gbar) h^mn dE_mn = target + total derivative",
                       zero_residual(qa.total_derivative_check(sp.simplify(r2["L2"] - tgt), s)), 0))
    recs.append(record("action.route1_vs_route2", "action",
                       "Routes 1 and 2 agree up to a total derivative",
                       zero_residual(qa.total_derivative_check(sp.simplify(r1["L2"] - r2["L2"]), s)), 0))
    recs.append(record("action.normalisation_control", "action",
                       "Negative control: doubling the target normalisation is detected",
                       zero_residual(qa.total_derivative_check(sp.simplify(r2["L2"] - 2 * tgt), s)), 0, expect_fail=True))
    el = qa.tensor_equation_from_action(s)
    eta, z = s["eta"], s["z"]
    fp = s["fp"]
    # E-L of (a^2/2)(f'^2 - f_z^2) is -(a^2 f')' + a^2 f_zz = -a^2 (f'' + 2 calH f' - f_zz)
    a2 = s["a"] ** 2
    calH = -1 / eta
    eom_target = -a2 * (sp.diff(fp, eta, 2) + 2 * calH * sp.diff(fp, eta) - sp.diff(fp, z, 2))
    recs.append(record("action.tensor_eom", "action",
                       "E-L equation of the tensor action is gamma'' + 2 calH gamma' + k^2 gamma = 0",
                       zero_residual(sp.simplify(el[0] - eom_target)), 0))
    fo = qa.first_order_variation_check()
    recs.append(record("action.first_order_total_derivative", "action",
                       "First-order bulk variation (generic diagonal h(eta)) is a total derivative on the background solution",
                       zero_residual(fo["correct"]["EL"]), 0))
    recs.append(record("action.first_order_wrong_lambda", "action",
                       "Negative control: with Lambda = 2H^2 the first-order variation is not a total derivative",
                       zero_residual(fo["wrong"]["EL"]), 0, expect_fail=True))
    return recs


def canonical_checks():
    recs = []
    z2, det = cv.derive_z_T_squared()
    eta = det["setup"]["eta"]
    H = det["setup"]["H"]
    M = det["M"]
    recs.append(record("canon.z_T_squared", "canonical",
                       "z_T^2 = M_pl^2 a^2/2 derived from the verified action", zero_residual(z2 - M**2 / (2 * H**2 * eta**2)), 0))
    recs.append(record("canon.effective_potential", "canonical", "z_T''/z_T = 2/eta^2",
                       zero_residual(cv.effective_potential() - 2 / sp.Symbol("eta", negative=True) ** 2), 0))
    el, _ = cv.verify_canonical_action()
    recs.append(record("canon.action_equivalence", "canonical",
                       "L[gamma = v/z] equals L_v up to a total derivative", zero_residual(el), 0))
    eq, e, v, zc = cv.mode_equation_from_action()
    k = sp.Symbol("k", positive=True)
    # Fourier substitute: v_zz -> -k^2 v
    expected = -(sp.diff(v, e, 2) - sp.diff(v, zc, 2) - 2 / e**2 * v)
    recs.append(record("canon.mode_equation", "canonical",
                       "E-L of the canonical action: v'' + (k^2 - 2/eta^2) v = 0", zero_residual(sp.simplify(eq - expected)), 0))
    C2, W = mn.normalization_from_wronskian()
    Hh, Mm = sp.symbols("H M_pl", positive=True)
    kk = sp.Symbol("k", positive=True)
    recs.append(record("canon.bd_normalisation", "canonical",
                       "Wronskian condition fixes |C|^2 = H^2/(M_pl^2 k^3)", zero_residual(C2 - Hh**2 / (Mm**2 * kk**3)), 0))
    recs.append(record("canon.exact_mode_solves_eom", "canonical",
                       "gamma = (1 + i k eta) e^{-i k eta} solves gamma'' - (2/eta) gamma' + k^2 gamma = 0",
                       zero_residual(mn.verify_mode_solves_equation()), 0))
    lim = mn.positive_frequency_limit()
    recs.append(record("canon.positive_frequency", "canonical",
                       "sqrt(2k) v_BD e^{i k eta} -> 1 as k eta -> -infinity", zero_residual(lim - 1), 0))
    return recs


def run_all():
    recs = []
    for f in (background_checks, metric_expansion_checks, linearized_einstein_checks, gauge_checks,
              svt_checks, action_checks, canonical_checks):
        recs += f()
    return recs
