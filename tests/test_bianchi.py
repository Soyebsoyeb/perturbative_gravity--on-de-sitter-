import sympy as sp
from conftest import select, assert_all_pass


def test_all_bianchi_records(symbolic_records):
    recs = select(symbolic_records, "bianchi.")
    assert_all_pass(recs)
    ids = {r.id for r in recs}
    assert {"bianchi.covariant_form", "bianchi.mixed_form", "bianchi.field_equation_divergence",
            "bianchi.naive_divergence_nonzero"} <= ids


def test_naive_divergence_is_a_negative_control(symbolic_records):
    r = [r for r in symbolic_records if r.id == "bianchi.naive_divergence_nonzero"][0]
    assert r.kind == "negative_control" and r.residual == 1.0 and r.passed
