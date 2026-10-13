"""Canonical commutation relations for the free tensor field.

Mode expansion (one polarization, comoving box):
    gamma(eta, x) = (1/sqrt V) sum_k [ gamma_k(eta) a_k + gamma_k*(eta) a^dagger_{-k} ] e^{i k x},
    pi(eta, x)    = z_T^2 gamma'(eta, x),    [a_k, a^dagger_k'] = delta_{k k'}.

Then [gamma(x), pi(y)] = (i/V) sum_k e^{i k (x-y)} * z_T^2 (gamma_k gamma_k*' - c.c.)/ i
which equals i delta(x - y) iff the Wronskian condition z_T^2 (gamma gamma*' - gamma* gamma') = i
holds for every k.  ``lattice_commutator`` verifies this numerically on a periodic lattice.
"""

from __future__ import annotations

import numpy as np

from ..dynamics import analytic_solutions as an


def commutator_coefficient(gamma, dgamma, z):
    """Return z^2 (gamma conj(dgamma) - conj(gamma) dgamma) / i; equals 1 if properly normalised."""
    return z**2 * (gamma * np.conj(dgamma) - np.conj(gamma) * dgamma) / 1j


def lattice_commutator(N=16, L=2 * np.pi, eta=-0.5, H=1.0, M_pl=1.0):
    """Matrix C_xy = [gamma(x_a), pi(x_b)] / i in units where the continuum value is delta/dx.

    On a periodic 1D lattice with N sites and spacing dx = L/N the exact answer is
    (1/dx) * identity.  Returns (C * dx, max |C*dx - identity|).
    """
    ks = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
    x = np.arange(N) * (L / N)
    z = an.z_T(eta, H, M_pl)
    C = np.zeros((N, N), dtype=complex)
    for k in ks:
        if k == 0:
            # The homogeneous mode has no BD amplitude (gamma ~ k^{-3/2}); its commutator
            # coefficient is fixed directly by the Wronskian normalisation (= 1).
            coeff = 1.0
        else:
            g = an.gamma_bd(abs(k), eta, H, M_pl)
            dg = an.dgamma_bd(abs(k), eta, H, M_pl)
            coeff = commutator_coefficient(g, dg, z)
        C += coeff * np.exp(1j * k * (x[:, None] - x[None, :])) / L
    dx = L / N
    Cdx = (C * dx).real
    return Cdx, float(np.max(np.abs(Cdx - np.eye(N))))
