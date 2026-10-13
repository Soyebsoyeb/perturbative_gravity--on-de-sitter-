"""Independent analytic benchmarks.

An analytic benchmark validates the *implementation* (equations, normalisations,
conventions are mutually consistent and reproduce known closed forms).  It cannot show
that the free-field de Sitter model is a complete description of gravity; see
docs/limitations.md.
"""

from __future__ import annotations

import numpy as np

from .records import record
from .. import background as bg
from ..dynamics import analytic_solutions as an
from ..spectra import polarization as pol
from ..spectra.tensor_power import tensor_power_from_mode, power_exact, power_superhorizon
from ..quantization import wronskian as wr
from ..quantization.commutation_relations import lattice_commutator


def _rel(a, b):
    return float(np.max(np.abs(a - b) / np.abs(b)))


def run():
    recs = []
    rng = np.random.default_rng(11)
    # scale factor and conformal-time relations
    t = np.linspace(-3, 8, 200)
    eta = bg.conformal_time_from_t(t)
    e1 = float(np.max(np.abs(bg.scale_factor_eta(eta) - bg.scale_factor_t(t)) / bg.scale_factor_t(t)))
    recs.append(record("bench.scale_factor_relation", "analytic", "a(eta(t)) = exp(H t)", e1, 1e-13, kind="benchmark",
                       relative_error=e1))
    e2 = float(np.max(np.abs(bg.t_from_conformal_time(eta) - t)))
    recs.append(record("bench.time_roundtrip", "analytic", "t -> eta -> t round trip", e2, 1e-12, kind="benchmark",
                       absolute_error=e2))
    # derivative relation: d eta/dt = 1/a (centred difference)
    tt = np.linspace(0, 4, 4001)
    de = np.gradient(bg.conformal_time_from_t(tt), tt)
    e3 = float(np.max(np.abs(de[2:-2] - 1.0 / bg.scale_factor_t(tt)[2:-2]) * bg.scale_factor_t(tt)[2:-2]))
    recs.append(record("bench.deta_dt", "analytic", "d eta/dt = 1/a (finite difference)", e3, 1e-5, kind="benchmark"))
    # exact mode solves the ODE (analytic second derivative)
    ks = np.geomspace(0.05, 50, 20)
    etas = -np.geomspace(1e-3, 1e3, 40)
    K, E = np.meshgrid(ks, etas)
    res = np.abs(an.residual_exact_conformal(K, E)) / (K**2 * np.abs(an.gamma_bd(K, E)))
    recs.append(record("bench.mode_solves_ode", "analytic", "Exact BD gamma satisfies gamma'' + 2 calH gamma' + k^2 gamma = 0",
                       float(res.max()), 1e-11, kind="benchmark", absolute_error=float(res.max()),
                       assumptions="analytic second derivative; relative to k^2|gamma|"))
    # canonical mode equation: v = z gamma solves v'' + (k^2 - 2/eta^2) v = 0 (analytic v and dv: check derivative consistency)
    h = 1e-4 / (K + 1.0 / np.abs(E))   # step scaled to the local time scale of the mode
    v = an.v_bd(K, E)
    dv_num = (an.v_bd(K, E + h) - an.v_bd(K, E - h)) / (2 * h)
    dv_an = an.dv_bd(K, E)
    scale = K * np.abs(v) + np.abs(v) / np.abs(E) + np.abs(dv_an)
    ev = float(np.max(np.abs(dv_num - dv_an) / scale))
    recs.append(record("bench.dv_consistency", "analytic", "Analytic dv_BD/d eta matches centred differences (step ~ 1e-4 of local time scale)", ev, 1e-6,
                       kind="benchmark"))
    # v = z gamma
    zT = an.z_T(E)
    rv = float(np.max(np.abs(np.abs(zT * an.gamma_bd(K, E)) - np.abs(an.v_bd(K, E))) / np.abs(an.v_bd(K, E))))
    recs.append(record("bench.v_equals_z_gamma", "analytic", "|z_T gamma_BD| = |v_BD| (same state up to a constant phase)",
                       rv, 1e-12, kind="benchmark", relative_error=rv))
    # asymptotics
    xs = np.array([1e2, 1e3, 1e4])
    for x in xs:
        etax = -x / 1.0
        g = an.gamma_bd(1.0, etax)
        gs = an.gamma_subhorizon(1.0, etax)
        rel = float(np.abs(abs(g) / abs(gs) - 1))
        recs.append(record(f"bench.subhorizon.x{int(x)}", "analytic",
                           f"Subhorizon asymptotics |gamma|/|gamma_WKB| -> 1 with error ~ 1/(2 x^2) (x = {int(x)})", rel,
                           1.0 / x**2, kind="benchmark", relative_error=rel,
                           assumptions="amplitude comparison; phase agrees up to the constant i"))
    for x in (1e-1, 1e-2, 1e-3):
        g = an.gamma_bd(1.0, -x)
        gs = an.gamma_superhorizon(1.0)
        rel = float(np.abs(abs(g) / abs(gs) - 1))
        recs.append(record(f"bench.superhorizon.x{x:g}", "analytic",
                           f"Superhorizon freeze-out |gamma| -> H/(M k^1.5) with error x^2/2 (x = {x:g})",
                           abs(rel - 0.5 * x**2), 0.5 * x**4 + 1e-14, kind="benchmark", relative_error=rel))
    # Wronskian
    et = -np.geomspace(1e-3, 1e2, 30)
    W = an.gamma_bd(1.0, et) * np.conj(an.dgamma_bd(1.0, et)) - np.conj(an.gamma_bd(1.0, et)) * an.dgamma_bd(1.0, et)
    ew = _rel(W, an.wronskian_gamma_exact(et))
    recs.append(record("bench.wronskian", "analytic", "gamma gamma*' - c.c. = 2i H^2 eta^2/M^2", ew, 1e-12, kind="benchmark",
                       relative_error=ew))
    ew2 = float(np.max(wr.bd_wronskian_deviation(1.3, et)))
    recs.append(record("bench.canonical_wronskian", "analytic", "z_T^2 (gamma gamma*' - c.c.) = i", ew2, 1e-12,
                       kind="benchmark", absolute_error=ew2))
    # power spectrum
    k = np.geomspace(0.1, 100, 15)
    P_mode = tensor_power_from_mode(k, an.gamma_bd(k, -0.5))
    P_ex = power_exact(k, -0.5)
    ep = _rel(P_mode, P_ex)
    recs.append(record("bench.power_from_mode", "analytic", "Delta^2 = 2 k^3|gamma|^2/pi^2 equals (2H^2/pi^2 M^2)(1+k^2 eta^2)",
                       ep, 1e-12, kind="benchmark", relative_error=ep))
    recs.append(record("bench.power_superhorizon", "analytic", "Delta^2 -> 2 H^2/(pi^2 M_pl^2) for k|eta| -> 0",
                       abs(float(power_exact(1.0, -1e-6)) / power_superhorizon() - 1), 1e-11, kind="benchmark"))
    # scaling with H and M
    ratio = power_superhorizon(2.0, 4.0) / power_superhorizon(1.0, 1.0)
    recs.append(record("bench.power_scaling", "analytic", "Delta^2 scales as H^2/M_pl^2",
                       abs(ratio - 4.0 / 16.0), 1e-14, kind="benchmark"))
    # polarization tensors
    worst = 0.0
    worst_c = 0.0
    worst_t = 0.0
    for _ in range(25):
        kvec = rng.standard_normal(3)
        worst = max(worst, pol.check_normalization(kvec))
        worst_c = max(worst_c, pol.check_completeness(kvec))
        t1, t2 = pol.check_transverse_traceless(kvec)
        worst_t = max(worst_t, t1, t2)
    recs.append(record("bench.polarization_norm", "analytic", "e^lam_ij e^lam'_ij = 2 delta (25 random directions)", worst, 1e-13, kind="benchmark"))
    recs.append(record("bench.polarization_completeness", "analytic", "sum_lam e_ij e_kl = 2 Lambda_{ij,kl}", worst_c, 1e-13, kind="benchmark"))
    recs.append(record("bench.polarization_tt", "analytic", "k_i e_ij = 0 and e_ii = 0", worst_t, 1e-13, kind="benchmark"))
    kv = [0.4, -1.1, 0.7]
    eR, eL = pol.circular_basis(kv)
    hR, hL = pol.helicity_eigenvalue(kv, eR), pol.helicity_eigenvalue(kv, eL)
    recs.append(record("bench.helicity", "analytic", "Circular polarizations have helicity +2 and -2",
                       max(abs(hR - 2), abs(hL + 2)), 1e-12, kind="benchmark"))
    # commutator
    cdx, err = lattice_commutator(N=16, eta=-0.7)
    recs.append(record("bench.lattice_commutator", "analytic",
                       "[gamma(x), pi(y)] = i delta_{xy}/dx on a periodic lattice with BD modes", err, 1e-12, kind="benchmark",
                       assumptions="zero mode normalised by the Wronskian condition"))
    return recs
