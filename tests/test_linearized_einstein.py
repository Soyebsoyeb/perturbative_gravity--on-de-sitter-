import sympy as sp
from desitter_perturbations.background import symbolic_background
from desitter_perturbations.symbolic.curvature import Geometry
from desitter_perturbations.symbolic.einstein_tensor import (linearized_einstein_covariant,
                                                            linearized_einstein_bruteforce)
from conftest import select, assert_all_pass


def test_covariant_matches_bruteforce_small_example():
    b = symbolic_background()
    eta, (x, y, z), H = b["eta"], b["spatial"], b["H"]
    geom = Geometry(b["coords_conformal"], b["g_conformal"])
    h = sp.Matrix([[x, 0, 0, 0], [0, eta * y, 0, 0], [0, 0, z, x], [0, 0, x, 1]])
    cov = linearized_einstein_covariant(geom, h, 3 * H**2)
    e0, e1 = linearized_einstein_bruteforce(b["coords_conformal"], b["g_conformal"], h, 3 * H**2)
    assert e0 == sp.zeros(4, 4)
    diff = (cov["dE"] - e1).applyfunc(lambda e: sp.cancel(sp.together(e)))
    assert diff == sp.zeros(4, 4)


def test_trace_part_dR():
    """dR = -h^{mn} Rbar_mn + gbar^{mn} dR_mn, with Rbar_mn = 3 H^2 gbar_mn gives dR = -3 H^2 h + gbar^{mn} dR_mn."""
    b = symbolic_background()
    eta, (x, y, z), H = b["eta"], b["spatial"], b["H"]
    geom = Geometry(b["coords_conformal"], b["g_conformal"])
    h = sp.Matrix([[x * eta, 0, 0, 0], [0, y, 0, 0], [0, 0, z, 0], [0, 0, 0, eta]])
    out = linearized_einstein_covariant(geom, h, 3 * H**2)
    trh = sum(geom.ginv[i, j] * h[i, j] for i in range(4) for j in range(4))
    lhs = out["dR"]
    rhs = -3 * H**2 * trh + sum(geom.ginv[i, j] * out["dRicci"][i, j] for i in range(4) for j in range(4))
    assert sp.simplify(lhs - rhs) == 0


def test_records(symbolic_records):
    assert_all_pass(select(symbolic_records, "lin."))
