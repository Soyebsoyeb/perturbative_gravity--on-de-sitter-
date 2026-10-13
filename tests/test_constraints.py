import sympy as sp
from desitter_perturbations.perturbations import constraints
from conftest import select, assert_all_pass


def test_scalar_constraint_chain():
    c = constraints.scalar_constraint_solution()
    k, ps = c["k"], c["ps"]
    funcs = list(c["e_traceless"].atoms(sp.core.function.AppliedUndef))
    assert len(funcs) == 2                                   # traceless equation involves Phi and Psi only
    other = [f for f in funcs if f != ps][0]
    assert sp.simplify(c["e_traceless"].subs(other, ps)) == 0   # and is proportional to their difference
    assert sp.simplify(c["phi_solution"] - ps) == 0           # Phi = Psi
    assert sp.simplify(c["e00_after_constraints"] + 2 * k**2 * ps) == 0   # Psi = 0 for k != 0


def test_constraint_records(symbolic_records):
    assert_all_pass(select(symbolic_records, "constraints."))
