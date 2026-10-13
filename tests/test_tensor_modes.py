import numpy as np
import pytest

from desitter_perturbations.dynamics import tensor_modes as tm, analytic_solutions as an
from desitter_perturbations.dynamics.numerical_solvers import SolverConfig, solve_mode
from desitter_perturbations.dynamics import residuals as rs
from desitter_perturbations.perturbations import tensor_sector as ts
from conftest import assert_all_pass


FORMS = ["conformal", "cosmic", "efold", "canonical"]


@pytest.mark.parametrize("form", FORMS)
def test_each_formulation_reproduces_exact_solution(form):
    s = solve_mode(1.0, 100.0, 1e-3, form, SolverConfig(rtol=1e-10, atol=1e-13, n_out=200))
    assert s.success
    assert np.max(rs.error_vs_exact(s)) < 5e-9


@pytest.mark.parametrize("k", [0.01, 1.0, 50.0])
def test_accuracy_independent_of_wavenumber(k):
    s = solve_mode(k, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-10, atol=1e-13, n_out=100))
    assert np.max(rs.error_vs_exact(s)) < 5e-9


def test_horizon_crossing_amplitude_is_sqrt2():
    s = solve_mode(1.0, 50.0, 1.0, "conformal", SolverConfig(rtol=1e-11, atol=1e-14, n_out=3))
    assert abs(abs(s.gamma[-1]) / an.gamma_ref(1.0) - np.sqrt(2.0)) < 1e-8


def test_subhorizon_decay_and_superhorizon_freeze():
    s = solve_mode(1.0, 200.0, 1e-3, "conformal", SolverConfig(rtol=1e-10, atol=1e-13, n_out=400))
    amp = np.abs(s.gamma)
    sub = s.x > 50
    assert np.allclose(amp[sub] * s.a[sub], amp[sub][0] * s.a[sub][0], rtol=5e-3)   # |gamma| ~ 1/a
    sup = s.x < 1e-2
    assert (amp[sup].max() - amp[sup].min()) / amp[sup].mean() < 1e-4              # frozen


def test_invalid_arguments():
    with pytest.raises(ValueError):
        solve_mode(1.0, 1.0, 10.0)
    with pytest.raises(ValueError):
        solve_mode(1.0, 10.0, 1.0, formulation="nonsense")


def test_rhs_conformal_matches_equation():
    k, eta = 2.0, -3.0
    y = tm.pack(0.3 + 0.1j, -0.2 + 0.5j)
    out = tm.rhs_conformal(eta, y, k)
    g, dg = 0.3 + 0.1j, -0.2 + 0.5j
    d2 = (2 / eta) * dg - k**2 * g
    assert np.allclose(out, [dg.real, dg.imag, d2.real, d2.imag])


def test_cosmic_and_conformal_equations_equivalent():
    """Transform a conformal-time state to cosmic time and check both RHS agree."""
    k, H, eta = 1.7, 1.0, -0.4
    a = -1 / (H * eta)
    t = np.log(a) / H
    g, dg = 0.2 - 0.3j, 0.5 + 0.1j            # gamma, gamma' (conformal)
    d2g = (2 / eta) * dg - k**2 * g
    gdot = dg / a
    # gamma_tt = (gamma'' - calH gamma') / a^2
    gtt_expected = (d2g - (-1 / eta) * dg) / a**2
    out = tm.rhs_cosmic(t, tm.pack(g, gdot), k, H)
    assert np.allclose(out[2] + 1j * out[3], gtt_expected)


def test_tensor_sector_symbolic_equation(symbolic_records):
    r = ts.tensor_dE()["residuals"]
    assert r["xx"] == 0 and r["yy"] == 0 and r["xy"] == 0 and all(v == 0 for v in r["others"])
