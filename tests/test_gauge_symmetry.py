import sympy as sp
from desitter_perturbations.symbolic import gauge_transformations as gt
from conftest import select, assert_all_pass


def test_svt_transformation_laws():
    assert all(v == 0 for v in gt.check_svt_laws().values())


def test_gauge_records(symbolic_records):
    assert_all_pass(select(symbolic_records, "gauge."))


def test_vector_gauge_function_does_not_touch_scalars():
    s = gt.svt_transformation_laws()
    # the transverse vector parameter Lx only appears in S_x and F_x shifts
    assert not s["dA"].has(s["Lx"]) and not s["dpsi"].has(s["Lx"])
    assert s["dS_x"].has(s["Lx"])
