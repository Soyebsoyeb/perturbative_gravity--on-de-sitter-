import pytest
import sympy as sp

from desitter_perturbations.symbolic import quadratic_action as qa
from desitter_perturbations.quantization import canonical_variables as cv
from conftest import select, assert_all_pass

pytestmark = pytest.mark.slow


def test_action_records(symbolic_records):
    assert_all_pass(select(symbolic_records, "action."))


def test_euler_expression_matches_known_case():
    t = sp.Symbol("t")
    f = sp.Function("f")(t)
    L = sp.Rational(1, 2) * sp.diff(f, t) ** 2 - sp.Rational(1, 2) * f**2     # harmonic oscillator
    assert sp.simplify(qa.euler_expression(L, f) + sp.diff(f, t, 2) + f) == 0
    # second-derivative Lagrangian: L = f''^2 gives f'''' (sign +)
    L2 = sp.diff(f, t, 2) ** 2
    assert sp.simplify(qa.euler_expression(L2, f) - 2 * sp.diff(f, t, 4)) == 0


def test_canonical_derivations(symbolic_records):
    assert_all_pass(select(symbolic_records, "canon."))
    z2, d = cv.derive_z_T_squared()
    assert sp.simplify(z2 - d["M"] ** 2 / (2 * d["setup"]["H"] ** 2 * d["setup"]["eta"] ** 2)) == 0
