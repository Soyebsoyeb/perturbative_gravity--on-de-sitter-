import pytest

from desitter_perturbations.validation import symbolic_checks


@pytest.fixture(scope="session")
def symbolic_records():
    """All symbolic verification records, computed once per test session."""
    return symbolic_checks.run_all()


def select(records, prefix):
    out = [r for r in records if r.id.startswith(prefix)]
    assert out, f"no records with prefix {prefix!r}"
    return out


def assert_all_pass(records):
    bad = [f"{r.id}: residual {r.residual:.3e} vs tol {r.tolerance:.1e}" for r in records if not r.passed]
    assert not bad, "failed checks:\n" + "\n".join(bad)
