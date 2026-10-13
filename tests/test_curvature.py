import sympy as sp
import pytest

from desitter_perturbations.background import symbolic_background
from desitter_perturbations.symbolic.curvature import Geometry
from desitter_perturbations.symbolic.einstein_tensor import einstein_tensor, field_equation


@pytest.fixture(scope="module")
def geoms():
    b = symbolic_background()
    return b, Geometry(b["coords_cosmic"], b["g_cosmic"]), Geometry(b["coords_conformal"], b["g_conformal"])


def test_ricci_scalar_and_invariants(geoms):
    b, gc, gf = geoms
    H = b["H"]
    for g in (gc, gf):
        assert sp.simplify(g.ricci_scalar - 12 * H**2) == 0
        assert sp.simplify(g.kretschmann() - 24 * H**4) == 0
        assert sp.simplify(g.ricci_squared() - 36 * H**4) == 0


def test_constant_curvature_identity(geoms):
    b, gc, gf = geoms
    for g in (gc, gf):
        assert all(v == 0 for v in g.constant_curvature_residual(b["H"]).values())


def test_einstein_and_field_equation(geoms):
    b, gc, gf = geoms
    H = b["H"]
    for g in (gc, gf):
        assert (einstein_tensor(g) + 3 * H**2 * g.g).applyfunc(sp.simplify) == sp.zeros(4, 4)
        assert field_equation(g, 3 * H**2) == sp.zeros(4, 4)


def test_wrong_lambda_is_not_a_solution(geoms):
    b, gc, gf = geoms
    E = field_equation(gf, 2 * b["H"] ** 2)
    assert E != sp.zeros(4, 4)


def test_riemann_sign_flip_detected(geoms):
    b, _, _ = geoms
    flipped = Geometry(b["coords_conformal"], b["g_conformal"], sign=-1)
    assert sp.simplify(flipped.ricci_scalar + 12 * b["H"] ** 2) == 0
    assert any(v != 0 for v in flipped.constant_curvature_residual(b["H"]).values())


def test_sphere_has_positive_curvature():
    th, ph = sp.symbols("theta phi", positive=True)
    g = Geometry((th, ph), sp.diag(1, sp.sin(th) ** 2))
    assert sp.simplify(g.ricci_scalar - 2) == 0   # unit 2-sphere: R = +2
