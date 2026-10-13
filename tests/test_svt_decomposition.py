import numpy as np
import pytest

from desitter_perturbations.symbolic import tensor_decomposition as td
from desitter_perturbations.perturbations import svt_numeric as sv, tensor_sector as ts
from desitter_perturbations.perturbations import gauge_invariants as gi
from conftest import select, assert_all_pass


def test_symbolic_identities(symbolic_records):
    assert_all_pass(select(symbolic_records, "svt."))
    assert_all_pass(select(symbolic_records, "sector."))


@pytest.mark.parametrize("N", [8, 12, 16])
def test_round_trip(N):
    assert sv.reconstruction_error(sv.random_symmetric_field(N, seed=N)) < 1e-13


def test_tt_part_is_transverse_and_traceless():
    g = sv.tt_part_field(16, seed=2)
    assert np.max(np.abs(sv.spectral_divergence(g))) < 1e-13
    assert np.max(np.abs(g[0, 0] + g[1, 1] + g[2, 2])) < 1e-13


def test_finite_difference_divergence_is_second_order():
    errs = []
    for N in (12, 24):
        g = sv.tt_part_field(N, seed=2, kmax=2)
        errs.append(np.max(np.abs(sv.finite_difference_divergence(g))))
    assert errs[1] < errs[0] / 3.0     # ~ 4x for second order


def test_zero_mode_handling():
    N = 8
    h = np.zeros((3, 3, N, N, N))
    h[:, :] = np.diag([1.0, 2.0, 3.0])[:, :, None, None, None]   # constant tensor
    comp, kv = sv.decompose_periodic(h)
    assert np.isclose(comp["psi"][0, 0, 0], -6.0 / 6.0 * N**3)  # FFT normalisation: sum = N^3 * mean
    assert abs(comp["E"][0, 0, 0]) == 0 and np.all(comp["F"][0, 0, 0] == 0)
    assert sv.reconstruction_error(h) < 1e-12


def test_tt_projector_properties():
    k = np.array([0.3, -1.2, 0.8])
    L = ts.tt_projector(k)
    assert np.allclose(np.einsum("ijkl,klmn->ijmn", L, L), L, atol=1e-13)    # idempotent
    A = np.random.default_rng(0).standard_normal((3, 3)); A = A + A.T
    g = ts.tt_project(A, k)
    assert np.allclose(g @ k, 0, atol=1e-13) and abs(np.trace(g)) < 1e-13


def test_projector_undefined_at_zero_k():
    with pytest.raises(ValueError):
        ts.transverse_projector([0, 0, 0])


def test_gauge_invariants_numeric():
    eta, res = gi.numeric_gauge_experiment()
    for n in ("shifted_1", "shifted_2"):
        assert np.max(np.abs(res[n]["Phi"] - res["original"]["Phi"])) < 1e-12
        assert np.max(np.abs(res[n]["Psi"] - res["original"]["Psi"])) < 1e-12
        assert np.max(np.abs(res[n]["A"] - res["original"]["A"])) > 1e-3


def test_covariant_tt_equation():
    r = ts.covariant_tensor_equation()
    assert all(v == 0 for v in r["divergence"]) and r["trace"] == 0
    assert all(v == 0 for v in r["equation"])
