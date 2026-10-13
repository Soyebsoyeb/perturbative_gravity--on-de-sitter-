import numpy as np
import pytest

from desitter_perturbations.dynamics import convergence as cv
from desitter_perturbations.dynamics.numerical_solvers import SolverConfig, solve_mode, solve_rk4_efold
from desitter_perturbations.dynamics import analytic_solutions as an, residuals as rs


def test_observed_order_helper():
    h = np.array([1.0, 0.5, 0.25])
    assert abs(cv.observed_order(h, 3 * h**2) - 2.0) < 1e-12


def test_rk4_is_fourth_order():
    h, e, p = cv.rk4_order_study()
    assert 3.9 < p < 4.1
    assert np.all(np.diff(e) < 0)


def test_error_decreases_with_tolerance():
    rows = cv.tolerance_study(rtols=(1e-5, 1e-7, 1e-9))
    errs = [r["max_err"] for r in rows]
    assert errs[0] > errs[1] > errs[2]
    assert all(r["max_err"] < 1e3 * r["rtol"] for r in rows)


@pytest.mark.parametrize("method", ["RK45", "DOP853", "Radau", "LSODA"])
def test_all_methods_converge(method):
    cfg = SolverConfig(method=method, rtol=1e-8, atol=1e-11, n_out=200)
    s = solve_mode(1.0, 60.0, 1e-3, "conformal", cfg)
    assert s.success and np.max(rs.error_vs_exact(s)) < 1e-4


def test_start_time_error_is_bogoliubov_mixing():
    rows = cv.start_x_study(x_starts=(10, 40, 160), x_end=1e-3)
    for r in rows:
        assert abs(r["ratio"] - r["predicted"]) < 1e-8
        assert abs(r["abs_beta"] - 1 / (2 * r["x_start"] ** 2)) < 1e-9
    devs = [abs(r["ratio"] - 1) for r in rows]
    assert devs[-1] < devs[0]


def test_no_error_spike_at_horizon_crossing():
    x, e = cv.error_trajectories(rtols=(1e-8,))[1e-8]
    assert e[(x > 0.3) & (x < 3)].max() < 2 * e[x < 0.05].max()


def test_residual_converges_with_output_density():
    peaks = []
    for n in (1000, 4000):
        s = solve_mode(1.0, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-12, atol=1e-15, n_out=n))
        peaks.append(rs.equation_residuals(s)[1][s.x < 50].max())
    assert peaks[1] < peaks[0] / 50          # 4^4 = 256 ideal


def test_cross_formulation_agreement():
    sols, errs = cv.formulation_comparison(rtol=1e-10)
    assert max(errs.values()) < 5e-9
    assert max(rs.cross_formulation(sols).values()) < 5e-9


@pytest.mark.slow
def test_mpmath_crosscheck():
    from desitter_perturbations.dynamics.numerical_solvers import solve_mode_mpmath
    _, g_mp = solve_mode_mpmath(1.0, 10.0, 1e-1, dps=30)
    s = solve_mode(1.0, 10.0, 1e-1, "conformal", SolverConfig(rtol=1e-12, atol=1e-15, n_out=3))
    assert abs(s.gamma[-1] - g_mp) / abs(g_mp) < 1e-9
