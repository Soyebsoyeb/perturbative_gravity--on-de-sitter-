import numpy as np
import pytest

from desitter_perturbations import background as bg
from desitter_perturbations.conventions import Units
from conftest import select, assert_all_pass


def test_scale_factor_conformal_relations():
    t = np.linspace(-2, 6, 50)
    eta = bg.conformal_time_from_t(t, H=1.7)
    assert np.allclose(bg.scale_factor_eta(eta, 1.7), bg.scale_factor_t(t, 1.7), rtol=1e-13)
    assert np.allclose(bg.t_from_conformal_time(eta, 1.7), t, atol=1e-12)
    assert np.allclose(bg.conformal_hubble(eta), bg.scale_factor_t(t, 1.7) * 1.7, rtol=1e-13)


def test_horizon_variable_and_wavelength():
    eta = np.array([-10.0, -1.0, -0.1])
    assert np.allclose(bg.x_of_eta(2.0, eta), 2.0 * np.abs(eta))
    assert np.allclose(bg.eta_of_x(2.0, bg.x_of_eta(2.0, eta)), eta)
    # at horizon crossing the physical wavelength is 2 pi / H
    k, H = 3.0, 0.8
    eta_c = bg.eta_of_x(k, 1.0)
    assert np.isclose(bg.physical_wavelength(k, eta_c, H), 2 * np.pi / H)


def test_eta_must_be_negative():
    with pytest.raises(ValueError):
        bg.t_from_conformal_time(0.5)


def test_units_derived_quantities():
    u = Units(H=2.0, M_pl=3.0)
    assert u.Lambda == pytest.approx(12.0)
    assert u.G == pytest.approx(1 / (8 * np.pi * 9.0))
    assert u.tensor_amplitude == pytest.approx(2 * 4.0 / (np.pi**2 * 9.0))


def test_symbolic_background_identities(symbolic_records):
    assert_all_pass(select(symbolic_records, "bg."))
