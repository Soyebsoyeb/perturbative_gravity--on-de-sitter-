import numpy as np
import pytest

from desitter_perturbations.spectra import tensor_power as tp, polarization as pol, parameter_sweeps as ps
from desitter_perturbations.spectra import mode_spectrum as msp
from desitter_perturbations.dynamics import analytic_solutions as an
from desitter_perturbations.perturbations.tensor_sector import tt_projector, polarization_basis


def test_superhorizon_amplitude():
    assert np.isclose(tp.power_superhorizon(1.0, 1.0), 2 / np.pi**2)
    assert np.isclose(tp.power_superhorizon(2.0, 0.5), 2 * 4 / (np.pi**2 * 0.25))


def test_power_from_mode_matches_exact_formula():
    k = np.geomspace(0.1, 100, 10)
    eta = -0.37
    assert np.allclose(tp.tensor_power_from_mode(k, an.gamma_bd(k, eta)), tp.power_exact(k, eta), rtol=1e-12)


def test_polarization_convention_controls_prefactor():
    """With a different polarization norm the same mode gives a proportionally different spectrum."""
    g = an.gamma_bd(1.0, -1e-3)
    p2 = tp.tensor_power_from_mode(1.0, g, n_pol=2, pol_norm=2.0)
    p1 = tp.tensor_power_from_mode(1.0, g, n_pol=2, pol_norm=1.0)
    assert np.isclose(p2 / p1, 2.0)


@pytest.mark.parametrize("seed", range(5))
def test_polarization_tensors(seed):
    k = np.random.default_rng(seed).standard_normal(3)
    assert pol.check_normalization(k) < 1e-13
    assert pol.check_completeness(k) < 1e-13
    t, tr = pol.check_transverse_traceless(k)
    assert t < 1e-13 and tr < 1e-13


def test_polarization_sum_reproduces_projector_trace():
    """sum_lam e_ij e_ij = 4 = 2 * 2 polarizations, the origin of the factor n_pol * N_e."""
    k = np.array([1.0, 2.0, -0.5])
    ep, ex = polarization_basis(k)
    assert np.isclose(np.sum(ep * ep) + np.sum(ex * ex), 4.0)
    L = tt_projector(k)
    assert np.isclose(np.einsum("ijij->", L), 2.0)       # trace of the TT projector = number of polarizations


def test_numerical_spectrum_matches_analytic():
    r = ps.sweep_wavenumber(n=7, k_min=0.2, k_max=60.0, eta_end=-1e-2, x_start=200.0)
    assert r["rel_err"].max() < 1e-8


def test_scale_invariance_superhorizon():
    r = ps.sweep_wavenumber(n=7, k_min=0.2, k_max=5.0, eta_end=-1e-3, x_start=200.0)
    spread = (r["numeric"].max() - r["numeric"].min()) / r["numeric"].mean()
    # the only k-dependence is the exact factor (1 + x_end^2); its spread is x_end,max^2 - x_end,min^2
    predicted = r["x_end"].max() ** 2 - r["x_end"].min() ** 2
    assert spread == pytest.approx(predicted, rel=1e-2)
    assert np.allclose(r["numeric"], tp.power_superhorizon() * (1 + r["x_end"] ** 2), rtol=1e-8)


def test_tilt_vanishes_superhorizon_and_blue_inside():
    k = np.geomspace(0.01, 300, 200)
    P = tp.power_exact(k, -0.01)
    nt = tp.tensor_tilt(k, P)
    assert abs(nt[k < 0.5].max()) < 1e-3
    assert nt[-1] > 1.0                    # Delta^2 ~ k^2 deep inside the horizon


def test_variance_integral():
    k = np.geomspace(1e-3, 1e1, 2000)
    P = np.full_like(k, 2.0)
    assert np.isclose(tp.variance_from_spectrum(k, P), 2.0 * np.log(k[-1] / k[0]), rtol=1e-6)


def test_horizon_crossing_amplitude():
    assert abs(msp.horizon_crossing_amplitude() - np.sqrt(2)) < 1e-6


def test_bogoliubov_power_formula():
    be = 0.2 * np.exp(1j * 0.9)
    al = np.sqrt(1 + abs(be) ** 2)
    assert np.isclose(tp.power_bogoliubov_superhorizon(al, be), abs(al - be) ** 2 * tp.power_superhorizon())
