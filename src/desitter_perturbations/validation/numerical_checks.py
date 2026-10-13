"""Numerical verification: solver accuracy, convergence, stability, consistency."""

from __future__ import annotations

import os
import numpy as np

from .records import record
from ..dynamics import analytic_solutions as an
from ..dynamics import convergence as cv
from ..dynamics import residuals as rs
from ..dynamics import initial_conditions as ic
from ..dynamics.numerical_solvers import SolverConfig, solve_mode, solve_mode_mpmath
from ..spectra import parameter_sweeps as ps
from ..spectra.tensor_power import tensor_power_from_mode, power_exact
from ..perturbations import svt_numeric as sv
from ..perturbations import gauge_invariants as gi
from ..quantization import wronskian as wr


def run(fast: bool = False):
    recs = []
    # 1. four formulations vs exact
    sols, errs = cv.formulation_comparison(k=1.0, rtol=1e-10)
    for f, e in errs.items():
        recs.append(record(f"num.exact_vs_{f}", "numerical", f"Numerical {f} solution vs exact BD (rtol 1e-10)",
                           e, 5e-9, kind="numerical", relative_error=e, method="DOP853"))
    cross = rs.cross_formulation(sols)
    for f, e in cross.items():
        recs.append(record(f"num.cross_{f}", "numerical", f"{f} vs conformal-time solution, pointwise", e, 5e-9,
                           kind="numerical", relative_error=e))
    # wavenumber independence of accuracy
    worst = 0.0
    for k in (0.01, 1.0, 100.0):
        s = solve_mode(k, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-10, atol=1e-13, n_out=200))
        worst = max(worst, float(np.max(rs.error_vs_exact(s))))
    recs.append(record("num.wavenumber_independence", "numerical",
                       "Accuracy is independent of k when expressed in x = k/(aH) (k = 0.01, 1, 100)", worst, 5e-9,
                       kind="numerical", relative_error=worst))
    # 2. tolerance scaling
    rows = cv.tolerance_study(rtols=(1e-6, 1e-8, 1e-10))
    bad = max(r["max_err"] / r["rtol"] for r in rows)
    recs.append(record("num.tolerance_scaling", "numerical", "Max error <= 1000 x rtol for rtol in {1e-6,1e-8,1e-10}", bad, 1000.0,
                       kind="numerical", notes="error/rtol ratio"))
    mono = all(rows[i]["max_err"] > rows[i + 1]["max_err"] for i in range(len(rows) - 1))
    recs.append(record("num.tolerance_monotone", "numerical", "Error decreases monotonically with rtol", 0.0 if mono else 1.0, 0.0, kind="numerical"))
    # 3. RK4 order
    h, e, p = cv.rk4_order_study()
    recs.append(record("num.rk4_order", "numerical", "Fixed-step RK4 observed order is 4", abs(p - 4.0), 0.1, kind="numerical",
                       notes=f"observed order {p:.3f}"))
    # 4. horizon crossing: no error spike
    traj = cv.error_trajectories(rtols=(1e-8,))
    x, e = traj[1e-8]
    near = e[(x > 0.3) & (x < 3.0)]
    far = e[x < 0.05]
    spike = float(near.max() / max(far.max(), 1e-300))
    recs.append(record("num.no_horizon_spike", "numerical",
                       "Error near horizon crossing (0.3 < x < 3) is no larger than the late-time plateau by more than 2x",
                       spike, 2.0, kind="numerical"))
    # 5. Bogoliubov prediction
    rows = cv.start_x_study(x_starts=(10, 50, 200), x_end=1e-3)
    worst = max(abs(r["ratio"] - r["predicted"]) for r in rows)
    recs.append(record("num.bogoliubov_prediction", "numerical",
                       "Power ratio of leading-WKB initial data equals the exact |alpha u + beta u*|^2 prediction",
                       worst, 1e-8, kind="numerical", absolute_error=worst))
    beta_err = max(abs(r["abs_beta"] - 1.0 / (2.0 * r["x_start"] ** 2)) for r in rows)
    recs.append(record("num.beta_scaling", "numerical", "|beta| = 1/(2 x_start^2) for leading-WKB initial data", beta_err, 1e-9,
                       kind="numerical"))
    # 6. spectrum
    r = ps.sweep_wavenumber(n=9, k_min=0.2, k_max=100.0, eta_end=-1e-2, x_start=300.0)
    recs.append(record("num.spectrum_vs_exact", "numerical", "Numerical Delta_t^2(k) at fixed eta vs (2H^2/pi^2M^2)(1+k^2 eta^2)",
                       float(r["rel_err"].max()), 1e-8, kind="numerical", relative_error=float(r["rel_err"].max())))
    r = ps.sweep_late_time()
    recs.append(record("num.late_time_freezeout", "numerical", "Delta_t^2 vs exact (1+x^2) for x_end in [1e-4, 3]",
                       float(r["rel_err"].max()), 1e-8, kind="numerical"))
    r = ps.sweep_H_Mpl()
    sc = float(np.max(np.abs(r["numeric"] / r["predicted"] - 1)))
    recs.append(record("num.H_M_scaling", "numerical", "Numerical power scales as H^2/M_pl^2 over a grid of (H, M_pl)", sc, 1e-8,
                       kind="numerical"))
    r = ps.sweep_bogoliubov()
    recs.append(record("num.bogoliubov_sweep", "numerical", "Excited-state power vs prediction for 13 phases of beta",
                       float(np.max(np.abs(r["numeric"] - r["predicted"]))), 1e-8, kind="numerical"))
    # 7. Wronskian conservation along the numerical solution
    s = solve_mode(1.0, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-11, atol=1e-14, n_out=500))
    wc = wr.conservation_from_solution(s)
    recs.append(record("num.wronskian_conservation", "numerical", "z^2 (gamma gamma*' - c.c.) = i along the numerical solution",
                       wc, 1e-8, kind="numerical", absolute_error=wc))
    # 8. residual diagnostic decays with output density (interpolation floor)
    peaks = []
    for n in (1000, 4000):
        sn = solve_mode(1.0, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-12, atol=1e-15, n_out=n))
        _, dyn = rs.equation_residuals(sn)
        peaks.append(float(np.max(dyn[sn.x < 50])))
    order = float(np.log(peaks[0] / peaks[1]) / np.log(4.0))
    recs.append(record("num.residual_order", "numerical", "Spline-derivative residual converges at ~4th order in output density",
                       abs(order - 4.0), 0.5, kind="numerical", notes=f"order {order:.2f}; peak residual {peaks[1]:.2e}"))
    # 9. mpmath high precision cross-check
    if not fast:
        xe, g_mp = solve_mode_mpmath(1.0, 10.0, 1e-1, dps=30)
        s2 = solve_mode(1.0, 10.0, 1e-1, "conformal", SolverConfig(rtol=1e-12, atol=1e-15, n_out=5))
        err = abs(s2.gamma[-1] - g_mp) / abs(g_mp)
        recs.append(record("num.mpmath_crosscheck", "numerical", "Double-precision DOP853 vs 30-digit Taylor integration (mpmath)", float(err),
                           1e-9, kind="numerical", relative_error=float(err)))
        err2 = abs(g_mp - an.gamma_bd(1.0, -1e-1)) / abs(g_mp)
        recs.append(record("num.mpmath_vs_exact", "numerical", "30-digit integration vs exact BD solution", float(err2), 1e-10,
                           kind="numerical", notes="tolerance limited by float conversion of the exact formula"))
    # 10. stability under initial x
    P = []
    for xs in (50.0, 100.0, 200.0, 400.0):
        s = solve_mode(1.0, xs, 1e-3, "conformal", SolverConfig(rtol=1e-10, atol=1e-13, n_out=5), ic_name="bd_exact")
        P.append(tensor_power_from_mode(1.0, s.gamma[-1]))
    spread = float((max(P) - min(P)) / np.mean(P))
    recs.append(record("num.stability_x_start", "numerical", "Exact-BD data give the same late-time power for x_start in {50..400}",
                       spread, 1e-8, kind="numerical"))
    # 11. SVT numeric
    err = max(sv.reconstruction_error(sv.random_symmetric_field(N, seed=s_)) for N, s_ in ((8, 1), (12, 2), (16, 3)))
    recs.append(record("num.svt_roundtrip", "numerical", "Spectral SVT decomposition/reconstruction round trip (N = 8, 12, 16)", err, 1e-13, kind="numerical"))
    g = sv.tt_part_field(16, seed=4)
    d1 = float(np.max(np.abs(sv.spectral_divergence(g))))
    tr = float(np.max(np.abs(g[0, 0] + g[1, 1] + g[2, 2])))
    recs.append(record("num.tt_transverse", "numerical", "Spectral TT part is transverse: d_j gamma_ij = 0", d1, 1e-13, kind="numerical"))
    recs.append(record("num.tt_traceless", "numerical", "Spectral TT part is traceless", tr, 1e-13, kind="numerical"))
    fd = float(np.max(np.abs(sv.finite_difference_divergence(g))))
    recs.append(record("num.tt_fd_not_transverse", "numerical",
                       "Finite-difference divergence of the spectral TT field is O(dx^2), not zero (role of the discretisation)",
                       fd, 1e-13, kind="numerical", expect_fail=True))
    kn, eE, eF = sv.noise_amplification()
    m = kn < 3.0
    sE = float(np.polyfit(np.log(kn[m]), np.log(eE[m]), 1)[0])
    sF = float(np.polyfit(np.log(kn[m]), np.log(eF[m]), 1)[0])
    recs.append(record("num.small_k_amplification_E", "numerical", "Error in E grows like k^-2 at small k", abs(sE + 2.0), 0.5,
                       kind="numerical", notes=f"slope {sE:.2f}"))
    recs.append(record("num.small_k_amplification_F", "numerical", "Error in F grows like k^-1 at small k", abs(sF + 1.0), 0.5,
                       kind="numerical", notes=f"slope {sF:.2f}"))
    # 12. gauge invariance numeric
    eta, res = gi.numeric_gauge_experiment()
    inv = max(float(np.max(np.abs(res[n_][q] - res["original"][q])))
              for n_ in ("shifted_1", "shifted_2") for q in ("Phi", "Psi"))
    dep = min(float(np.max(np.abs(res[n_][q] - res["original"][q])))
              for n_ in ("shifted_1", "shifted_2") for q in ("A", "psi", "B", "E"))
    recs.append(record("num.gauge_invariants", "numerical", "Bardeen potentials are unchanged by gauge shifts", inv, 1e-12, kind="numerical",
                       absolute_error=inv))
    recs.append(record("num.gauge_dependents", "numerical", "A, psi, B, E do change under the same shifts (control)", dep, 1e-6,
                       kind="numerical", expect_fail=True))
    return recs
