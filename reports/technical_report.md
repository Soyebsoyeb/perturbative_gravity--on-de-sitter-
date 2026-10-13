# Technical Report: Perturbative Gravity on de Sitter Spacetime

## 1. Scope

This report summarizes the derivations, the verification evidence and the physical results of the `desitter_perturbations` package. It distinguishes established results, independently verified results, exploratory observations and new results. No new physics is claimed and no claim is made about solving quantum gravity.

## 2. Background and conventions

Signature $(-+++)$, $R^\rho{}_{\sigma\mu\nu}=\partial_\mu\Gamma^\rho_{\nu\sigma}-\dots$, $R_{\sigma\nu}=R^\rho{}_{\sigma\rho\nu}$. Flat-slicing de Sitter, $a=-1/(H\eta)$, satisfies

$$
R_{\mu\nu\rho\sigma}=H^2\big(g_{\mu\rho}g_{\nu\sigma}-g_{\mu\sigma}g_{\nu\rho}\big),\qquad R=12H^2,\qquad G_{\mu\nu}=-3H^2g_{\mu\nu},\qquad \Lambda=3H^2 .
$$

This is verified in cosmic and conformal slicings on all 256 Riemann components, with Kretschmann scalar $24H^4$ and Ricci-squared $36H^4$, and with a flipped-sign negative control.

## 3. Metric perturbation

$g=\bar g+\epsilon h$. The inverse metric and $\sqrt{-g}$ are expanded to fourth order by two methods; the measured truncation slopes are $0.99,\,2.01,\,3.01$ (inverse metric) and $0.97,\,2.00,\,3.01$ ($\sqrt{-g}$). The covariant $\delta\Gamma$ agrees with direct variation of the Christoffel symbols, and the Palatini identity holds.

## 4. Linearized equations and Bianchi identities

The linearized Einstein tensor from the covariant formula agrees with a brute-force expansion of the full nonlinear tensor. The covariant, mixed-index and divergence forms of the linearized Bianchi identity are verified. The naive coordinate divergence of $\delta E_{\mu\nu}$ is nonzero, which confirms the need for the covariant formulation.

## 5. Gauge structure

Pure-gauge perturbations $h=-\mathcal L_\xi\bar g$ solve the linearized equations identically. The SVT laws

$$
A\to A-(\mathcal HT+T'),\quad B\to B+T-L',\quad \psi\to\psi+\mathcal HT,\quad E\to E-L,\quad S_i\to S_i-L_i',\quad F_i\to F_i-L_i
$$

were derived by component matching and checked against the Lie derivative. Residual freedom: none in Newtonian gauge, $T=C/a$ in synchronous gauge, none in spatially flat gauge. Degree-of-freedom count: $10-4-4=2$.

## 6. Sector equations

Scalar: $\delta E_{00}$, $\delta E_{0z}$ and $\delta E_{xx}-\delta E_{zz}$ expressed through the Bardeen potentials. Vector: $-\tfrac12\partial_z^2V$ and $-\tfrac12\partial_z\big(V'+2\mathcal HV\big)$. Tensor:

$$
\delta E_{ij}=\tfrac12\big(\gamma_{ij}''+2\mathcal H\gamma_{ij}'-\nabla^2\gamma_{ij}\big),
$$

and in covariant TT form $\delta E_{\mu\nu}=-\tfrac12\big(\Box h_{\mu\nu}-2H^2h_{\mu\nu}\big)$.

## 7. Second-order action

Route one expands the exact action to second order; route two evaluates $-\tfrac12\sqrt{-\bar g}\,h^{\mu\nu}\delta E_{\mu\nu}$. They agree up to a total derivative, tested with a custom Euler-Lagrange operator after the off-the-shelf routine proved vacuous. The result is

$$
S_2=\frac{M_{\rm pl}^2}{8}\int d\eta\,d^3x\,a^2\big[(\gamma_{ij}')^2-(\partial_k\gamma_{ij})^2\big].
$$

## 8. Canonical quantization

$$
z_T^2=\frac{M_{\rm pl}^2a^2}{2}\ (e_{ij}e_{ij}=2),\qquad \frac{z_T''}{z_T}=\frac{2}{\eta^2},\qquad \pi=z_T^2\gamma',\qquad z_T^2\big(\gamma\gamma^{*\prime}-\gamma^*\gamma'\big)=i,
$$

with Bunch-Davies mode $v_k=\dfrac{e^{-ik\eta}}{\sqrt{2k}}\Big(1-\dfrac{i}{k\eta}\Big)$.

## 9. Spectrum

$$
\Delta_t^2=\frac{2H^2}{\pi^2M_{\rm pl}^2}\,(1+x^2),\qquad x=\frac{k}{aH}.
$$

Numerical agreement: $9.9\times10^{-12}$.

## 10. Numerical evidence

| Check | Value |
|---|---|
| Four formulations vs exact, $\mathrm{rtol}=10^{-10}$ | $5.5\times10^{-10}$ to $7.9\times10^{-10}$ |
| Cross-formulation | about $3.1\times10^{-10}$ |
| RK4 order | $4.003$ |
| mpmath cross-check | $2.6\times10^{-13}$ |
| Wronskian drift | $6.3\times10^{-11}$ |
| Bogoliubov prediction | $1.9\times10^{-11}$ |
| Gauge-invariant change | $3.6\times10^{-15}$ |
| Residual convergence order | $4.00$ |

## 11. Exploratory observation

Leading-WKB initial data at finite $x_{\rm start}$ is a Bogoliubov state with $\lvert\beta\rvert=1/(2x_{\rm start}^2)$; the power ratio follows the exact prediction. This is an exact consequence of the mode equation and is reported as a numerical-practice note.

## 12. Verification summary

128 of 128 records pass, including 9 negative controls, on Python 3.13.16, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0. See [validation_report.md](validation_report.md) and `validation_report.json`.

## 13. Limitations

Free field on exact de Sitter; no interactions, loops or matter; Bunch-Davies is an input; the polarization normalization affects numerical factors. See [../docs/limitations.md](../docs/limitations.md).
