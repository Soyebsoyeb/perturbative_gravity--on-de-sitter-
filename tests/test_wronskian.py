import numpy as np

from desitter_perturbations.quantization import wronskian as wr, commutation_relations as cr, vacuum_states as vs
from desitter_perturbations.quantization import mode_normalization as mn
from desitter_perturbations.dynamics import analytic_solutions as an
from desitter_perturbations.dynamics.numerical_solvers import SolverConfig, solve_mode


def test_canonical_wronskian_of_bd_mode_is_i():
    k = 0.8
    eta = np.array([-30.0, -3.0, -0.1])
    W = wr.canonical_wronskian(an.v_bd(k, eta), an.dv_bd(k, eta))
    assert np.allclose(W, 1j, atol=1e-12)


def test_gamma_wronskian_bd():
    eta = -np.geomspace(1e-3, 1e2, 20)
    assert np.max(wr.bd_wronskian_deviation(1.1, eta)) < 1e-12


def test_wronskian_conserved_by_numerical_solution():
    s = solve_mode(1.0, 100.0, 1e-3, "conformal", SolverConfig(rtol=1e-11, atol=1e-14, n_out=300))
    assert wr.conservation_from_solution(s) < 1e-8


def test_bogoliubov_preserves_wronskian():
    k, eta = 1.0, -2.0
    be = 0.4 * np.exp(0.7j)
    al = np.sqrt(1 + abs(be) ** 2)
    g = vs.bogoliubov_mode(k, eta, al, be)
    ug, dug = an.gamma_bd(k, eta), an.dgamma_bd(k, eta)
    dg = al * dug + be * np.conj(dug)
    assert abs(wr.gamma_wronskian(g, dg, an.z_T(eta)) - 1j) < 1e-12


def test_unnormalised_mode_violates_commutator():
    k, eta = 1.0, -2.0
    g, dg = 2.0 * an.gamma_bd(k, eta), 2.0 * an.dgamma_bd(k, eta)    # wrong normalisation
    assert abs(cr.commutator_coefficient(g, dg, an.z_T(eta)) - 1) > 1.0


def test_lattice_commutator_is_identity():
    _, err = cr.lattice_commutator(N=12, eta=-0.9)
    assert err < 1e-12


def test_symbolic_normalisation():
    import sympy as sp
    C2, _ = mn.normalization_from_wronskian()
    H, M, k = sp.symbols("H M_pl k", positive=True)
    assert sp.simplify(C2 - H**2 / (M**2 * k**3)) == 0
    assert mn.positive_frequency_limit() == 1


def test_vacuum_state_helpers():
    assert vs.check_bogoliubov(np.sqrt(1.09), 0.3)
    assert not vs.check_bogoliubov(1.0, 0.3)
    assert abs(vs.wkb_state_beta(10.0) - 0.005) < 1e-15
    assert np.isclose(vs.superhorizon_power_factor(1.0, 0.0), 1.0)
