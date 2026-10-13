"""Conventions used throughout the package.

Every convention-dependent statement in the code and documentation refers to
the choices recorded here.  They are checked against the implementation by
``validation.symbolic_checks`` rather than assumed.

Summary
-------
* Signature (-,+,+,+), four spacetime dimensions, units c = hbar = 1.
* Christoffel symbols  Gamma^l_{mn} = (1/2) g^{lr}(d_m g_{rn} + d_n g_{rm} - d_r g_{mn}).
* Riemann tensor       R^r_{s m n} = d_m Gamma^r_{n s} - d_n Gamma^r_{m s}
                                     + Gamma^r_{m l} Gamma^l_{n s} - Gamma^r_{n l} Gamma^l_{m s}.
* Ricci tensor         R_{s n} = R^r_{s r n};  Ricci scalar R = g^{s n} R_{s n}.
  With these choices a sphere has positive curvature and de Sitter space has
  R_{mnrs} = H^2 (g_{mr} g_{ns} - g_{ms} g_{nr}),  R = 12 H^2.
* Field equation       E_{mn} = G_{mn} + Lambda g_{mn} = 0,  Lambda = 3 H^2.
* Action               S = (1/(16 pi G)) int d^4x sqrt(-g) (R - 2 Lambda),
  reduced Planck mass  M_pl^2 = 1/(8 pi G).
* Perturbation         g_{mn} = gbar_{mn} + eps h_{mn}; indices of h are raised
  with the background metric.
* Gauge transformation h_{mn} -> h_{mn} - nabla_m xi_n - nabla_n xi_m  (nabla is
  the background connection).  This is the passive transformation
  h -> h - L_xi gbar, i.e. the change of the perturbation under the
  coordinate change x'^m = x^m + eps xi^m (see docs/gauge_symmetry.md).
* Fourier transform    f(x) = int d^3k/(2 pi)^3 f_k exp(i k.x).
* Tensor perturbation  ds^2 = a^2 [ -d eta^2 + (delta_ij + gamma_ij) dx^i dx^j ].
* Polarization tensors e^lam_ij e^lam'_ij = 2 delta^{lam lam'},  gamma_ij = sum_lam e^lam_ij gamma_lam.
* Tensor spectrum      <gamma_ij(x) gamma_ij(x)> = int dln k Delta_t^2(k).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math

SIGNATURE = (-1, 1, 1, 1)
DIMENSION = 4
RIEMANN_SIGN = +1  # +1 selects the convention documented above


@dataclass(frozen=True)
class Units:
    """Numerical unit system.  Everything is dimensionless once H and M_pl are set.

    ``H`` is the Hubble rate and ``M_pl`` the reduced Planck mass, both in the
    same (arbitrary) energy unit.  Wavenumbers are measured in units of H.
    """

    H: float = 1.0
    M_pl: float = 1.0

    @property
    def Lambda(self) -> float:
        return 3.0 * self.H**2

    @property
    def G(self) -> float:
        return 1.0 / (8.0 * math.pi * self.M_pl**2)

    @property
    def tensor_amplitude(self) -> float:
        """Exact-de-Sitter superhorizon tensor spectrum 2 H^2 / (pi^2 M_pl^2)."""
        return 2.0 * self.H**2 / (math.pi**2 * self.M_pl**2)

    def as_dict(self) -> dict:
        d = asdict(self)
        d.update(Lambda=self.Lambda, G=self.G)
        return d


DEFAULT_UNITS = Units()

# Polarization normalisation: e^lam_ij e^lam'_ij = POLARIZATION_NORM * delta^{lam lam'}
POLARIZATION_NORM = 2.0

CONVENTIONS_TABLE = {
    "signature": "(-,+,+,+)",
    "dimension": 4,
    "units": "c = hbar = 1",
    "riemann": "R^r_{s m n} = d_m Gamma^r_{n s} - d_n Gamma^r_{m s} + Gamma^r_{m l} Gamma^l_{n s} - Gamma^r_{n l} Gamma^l_{m s}",
    "ricci": "R_{s n} = R^r_{s r n}",
    "einstein_equation": "G_{mn} + Lambda g_{mn} = 0, Lambda = 3 H^2",
    "action": "S = (1/16 pi G) int sqrt(-g) (R - 2 Lambda)",
    "reduced_planck_mass": "M_pl^2 = 1/(8 pi G)",
    "perturbation": "g = gbar + eps h",
    "gauge": "h -> h - nabla_m xi_n - nabla_n xi_m",
    "fourier": "f(x) = int d^3k/(2 pi)^3 f_k exp(i k.x)",
    "polarization_norm": "e^lam_ij e^lam'_ij = 2 delta",
    "tensor_spectrum": "<gamma_ij gamma_ij> = int dlnk Delta_t^2",
}
