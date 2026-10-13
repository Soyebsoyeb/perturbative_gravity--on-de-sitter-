# Conventions

Every convention-dependent statement in the code refers to this page. The conventions are not assumed: each one is exercised by a verification record (identifiers in brackets, see `reports/validation_report.md`).

## Spacetime and units

* Four spacetime dimensions, signature $(-,+,+,+)$, units $c=\hbar=1$.
* Newton's constant $G$ and reduced Planck mass $M_{\rm pl}^2=1/(8\pi G)$.
* Numerical work uses units with $H$ and $M_{\rm pl}$ as free parameters (`conventions.Units`). Wavenumbers are measured in the same units as $H$ with $a(t=0)=1$ in the cosmic-time formulations. Results are expressed through the dimensionless horizon variable $x=k/(aH)$ wherever possible, so they do not depend on this choice.

## Curvature

$$
\Gamma^\lambda{}_{\mu\nu}=\tfrac12 g^{\lambda\rho}(\partial_\mu g_{\rho\nu}+\partial_\nu g_{\rho\mu}-\partial_\rho g_{\mu\nu}),
$$
$$
R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho{}_{\nu\sigma}-\partial_\nu\Gamma^\rho{}_{\mu\sigma}+\Gamma^\rho{}_{\mu\lambda}\Gamma^\lambda{}_{\nu\sigma}-\Gamma^\rho{}_{\nu\lambda}\Gamma^\lambda{}_{\mu\sigma},
\qquad R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu},\qquad R=g^{\sigma\nu}R_{\sigma\nu}.
$$

With these choices the unit two-sphere has $R=+2$ (test `test_sphere_has_positive_curvature`) and de Sitter space satisfies

$$
\bar R_{\mu\nu\rho\sigma}=H^2(\bar g_{\mu\rho}\bar g_{\nu\sigma}-\bar g_{\mu\sigma}\bar g_{\nu\rho}),\quad
\bar R_{\mu\nu}=3H^2\bar g_{\mu\nu},\quad \bar R=12H^2,\quad \bar G_{\mu\nu}=-3H^2\bar g_{\mu\nu}.
$$

[`bg.constant_curvature.*`, `bg.ricci_scalar.*`, `bg.einstein.*`]. A deliberately flipped Riemann sign gives $R=-12H^2$ and is flagged by the negative control [`bg.sign_control`].

## Field equation and action

$$
E_{\mu\nu}\equiv G_{\mu\nu}+\Lambda g_{\mu\nu}=0,\qquad \Lambda=3H^2,\qquad
S=\frac{1}{16\pi G}\int d^4x\,\sqrt{-g}\,(R-2\Lambda).
$$

Variation gives $\delta S=-\frac{1}{16\pi G}\int d^4x\sqrt{-g}\,E^{\mu\nu}\delta g_{\mu\nu}$ up to boundary terms. A wrong $\Lambda$ is detected by the first-order variation check [`action.first_order_wrong_lambda`].

## Perturbation

$$
g_{\mu\nu}=\bar g_{\mu\nu}+\epsilon\,h_{\mu\nu}.
$$

The bookkeeping parameter $\epsilon$ is kept until orders have been identified. Indices of $h$ are raised and lowered with $\bar g$: $h^{\mu\nu}=\bar g^{\mu\alpha}\bar g^{\nu\beta}h_{\alpha\beta}$, $h=\bar g^{\mu\nu}h_{\mu\nu}$. Expansions of the inverse metric and $\sqrt{-g}$ are in `docs/metric_variation.md`.

## Gauge transformation

$$
h_{\mu\nu}\ \to\ h_{\mu\nu}-\bar\nabla_\mu\xi_\nu-\bar\nabla_\nu\xi_\mu = h_{\mu\nu}-(\mathcal L_\xi\bar g)_{\mu\nu},\qquad \xi_\mu=\bar g_{\mu\nu}\xi^\nu .
$$

This is the passive transformation induced by the coordinate change $x'^\mu=x^\mu+\epsilon\xi^\mu$. See `docs/gauge_symmetry.md` for the sign discussion [`gauge.lie_equals_covariant`].

## Coordinates for flat de Sitter slicing

| | cosmic time | conformal time |
|---|---|---|
| metric | $-dt^2+a^2d\mathbf x^2$ | $a^2(-d\eta^2+d\mathbf x^2)$ |
| scale factor | $a=e^{Ht}$ | $a=-1/(H\eta)$, $\eta<0$ |
| relation | $\eta=-e^{-Ht}/H$ | $d\eta=dt/a$ |
| Hubble rate | $H=\dot a/a$ | $\mathcal H=a'/a=-1/\eta=aH$ |

## Fourier transform and polarizations

$$
f(\mathbf x)=\int\frac{d^3k}{(2\pi)^3}f_{\mathbf k}\,e^{i\mathbf k\cdot\mathbf x},\qquad \partial_j\to ik_j .
$$

Tensor perturbation: $ds^2=a^2[-d\eta^2+(\delta_{ij}+\gamma_{ij})dx^idx^j]$, $\gamma_{ij}=\sum_\lambda e^\lambda_{ij}(\hat k)\gamma_\lambda(\mathbf k)$ with real linear polarizations normalised as

$$
e^\lambda_{ij}e^{\lambda'}_{ij}=2\delta^{\lambda\lambda'},\quad \hat k_ie^\lambda_{ij}=0,\quad e^\lambda_{ii}=0,\quad \sum_\lambda e^\lambda_{ij}e^\lambda_{kl}=2\Lambda_{ij,kl}.
$$

[`bench.polarization_*`]. The factor 2 is not cosmetic: it is what makes $z_T^2=M_{\rm pl}^2a^2/2$ and fixes the prefactor of the power spectrum.

## Tensor power spectrum

$$
\langle\gamma_{ij}(\mathbf x)\gamma_{ij}(\mathbf x)\rangle=\int d\ln k\,\Delta_t^2(k),\qquad
\Delta_t^2=\frac{n_{\rm pol}N_e\,k^3}{2\pi^2}|\gamma_k|^2=\frac{2k^3}{\pi^2}|\gamma_k|^2\quad(n_{\rm pol}=N_e=2).
$$

## Boundary terms

Total derivatives are discarded in the quadratic action (spatial periodicity or compact support; no assumption is made about the time boundaries beyond the existence of a well-posed variational problem). Every integration by parts is listed in `docs/quadratic_action.md`.
