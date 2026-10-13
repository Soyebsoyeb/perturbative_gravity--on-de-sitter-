import numpy as np
import pytest

from desitter_perturbations.dynamics import analytic_solutions as an, initial_conditions as ic
from desitter_perturbations.validation import analytic_benchmarks
from conftest import assert_all_pass


def test_analytic_benchmark_records():
    assert_all_pass(analytic_benchmarks.run())


def test_exact_solution_satisfies_equation():
    k = np.array([0.1, 1.0, 10.0])[:, None]
    eta = -np.geomspace(1e-2, 1e2, 30)[None, :]
    r = an.residual_exact_conformal(k, eta) / (k**2 * np.abs(an.gamma_bd(k, eta)))
    assert np.max(np.abs(r)) < 1e-11


def test_finite_difference_derivative_matches_analytic():
    k, eta, h = 1.3, -2.0, 1e-5
    fd = (an.gamma_bd(k, eta + h) - an.gamma_bd(k, eta - h)) / (2 * h)
    assert abs(fd - an.dgamma_bd(k, eta)) < 1e-8


def test_amplitude_formula():
    k, eta = 0.7, -0.3
    assert np.isclose(abs(an.gamma_bd(k, eta)) ** 2, an.abs_gamma_squared_exact(k, eta))


def test_scaling_with_H_and_M():
    a1 = abs(an.gamma_bd(1.0, -1e-6, H=1.0, M_pl=1.0))
    a2 = abs(an.gamma_bd(1.0, -1e-6, H=3.0, M_pl=2.0))
    assert np.isclose(a2 / a1, 3.0 / 2.0, rtol=1e-9)


def test_bogoliubov_roundtrip():
    k, eta0 = 1.0, -30.0
    al, be = 1.2, 0.3 + 0.4j
    al = np.sqrt(1 + abs(be) ** 2)
    g, dg = ic.bogoliubov(k, eta0, al, be)
    al2, be2 = ic.bogoliubov_coefficients(k, eta0, g, dg)
    assert abs(al2 - al) < 1e-12 and abs(be2 - be) < 1e-12


def test_bd_coefficients_trivial():
    al, be = ic.bogoliubov_coefficients(1.0, -10.0, *ic.bd_exact(1.0, -10.0))
    assert abs(al - 1) < 1e-13 and abs(be) < 1e-13


def test_wkb_data_has_beta_one_over_2x_squared():
    for x0 in (10.0, 50.0):
        al, be = ic.bogoliubov_coefficients(1.0, -x0, *ic.bd_asymptotic(1.0, -x0))
        assert abs(abs(be) - 1 / (2 * x0**2)) < 1e-12
        assert abs(abs(al) ** 2 - abs(be) ** 2 - 1) < 1e-12
