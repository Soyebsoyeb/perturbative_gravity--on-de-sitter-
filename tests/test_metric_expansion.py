import sympy as sp
from desitter_perturbations.symbolic import metric_expansion as me
from conftest import select, assert_all_pass


def test_inverse_metric_second_order_diagonal_generic():
    res = me.check_inverse_metric(me.generic_diagonal_background(), me.generic_perturbation(), 2)
    assert all(e == 0 for e in res)


def test_inverse_metric_mutations_detected():
    gbar, h = me.generic_diagonal_background(), me.generic_perturbation()
    bad = me.check_inverse_metric(gbar, h, 2, builder=me.inverse_metric_series_wrong_sign)
    assert any(e != 0 for e in bad)


def test_sqrt_det_formula_and_mutation():
    gbar, h = me.generic_diagonal_background(), me.generic_perturbation()
    assert me.check_sqrt_det(gbar, h) == (0, 0)
    d1, d2 = me.check_sqrt_det(gbar, h, formula=me.sqrt_det_series_formula_wrong_factor)
    assert d1 == 0 and d2 != 0


def test_two_dimensional_example():
    g, h = me.two_dimensional_example()
    assert all(e == 0 for e in me.check_inverse_metric(g, h, 2))
    assert me.check_sqrt_det(g, h) == (0, 0)


def test_truncate_orders():
    e = me.EPS
    assert me.truncate(1 + e + e**2 + e**3, 2) == 1 + e + e**2


def test_records(symbolic_records):
    assert_all_pass(select(symbolic_records, "exp."))
