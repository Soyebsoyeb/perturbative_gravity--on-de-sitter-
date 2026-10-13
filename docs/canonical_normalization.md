# Canonical Normalization of the Tensor Sector

Status: derived symbolically (`quantization.canonical_variables`) and verified by the `quant.*` records.

## Starting point

The exact second-order tensor action on de Sitter, derived by two independent routes in [quadratic_action.md](quadratic_action.md), is

$$
S_2=\frac{M_{\rm pl}^2}{8}\int d\eta\,d^3x\;a^2\Big[(\gamma_{ij}')^2-(\partial_k\gamma_{ij})^2\Big],
$$

with $\gamma_{ij}$ transverse-traceless and $M_{\rm pl}^2=1/(8\pi G)$.

## Decomposition into polarizations

With $\gamma_{ij}=\sum_\lambda e^\lambda_{ij}\gamma_\lambda$ and the normalization $e^\lambda_{ij}e^{\lambda'}_{ij}=2\delta^{\lambda\lambda'}$, each polarization has

$$
S_2^{(\lambda)}=\frac{M_{\rm pl}^2}{4}\int d\eta\,d^3x\;a^2\Big[(\gamma_\lambda')^2-(\partial\gamma_\lambda)^2\Big].
$$

## Canonical variable

Define $v=z_T\gamma_\lambda$. Matching to the canonical form $\tfrac12\big[v'^2-(\nabla v)^2+\tfrac{z_T''}{z_T}v^2\big]$ up to a total derivative gives

$$
z_T^2=\frac{M_{\rm pl}^2a^2}{2},\qquad z_T=\frac{M_{\rm pl}\,a}{\sqrt2}.
$$

The code derives this by matching coefficients (`derive_z_T_squared`) rather than assuming it, and `verify_canonical_action` checks that the action in terms of $v$ equals the original up to a total derivative.

## Effective potential

On de Sitter $a=-1/(H\eta)$, so $z_T''/z_T=a''/a=2/\eta^2$ and the Mukhanov-Sasaki-type equation is

$$
v''+\Big(k^2-\frac{2}{\eta^2}\Big)v=0 .
$$

This is a symbolic identity checked in `quant.mode_equation_from_action`.

## Conjugate momentum

Per polarization, $\pi=\partial\mathcal L/\partial\gamma'=z_T^2\gamma'$. The Poisson bracket $\{\gamma,\pi\}=\delta$ is the starting point for the commutation relations in [quantum_state.md](quantum_state.md).

## Normalization convention caveat

The factor $2$ in $z_T^2$ depends entirely on the polarization normalization $e_{ij}e_{ij}=2$. With $e_{ij}e_{ij}=1$ one would obtain $z_T^2=M_{\rm pl}^2a^2$ and the power spectrum would differ by a factor of $2$ per polarization. The constant is stored in `conventions.POLARIZATION_NORM` and used everywhere.
