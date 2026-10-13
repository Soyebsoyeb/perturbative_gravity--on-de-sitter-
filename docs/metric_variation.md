# Metric variation

## Inverse metric

Write $g=\bar g(1+\epsilon M)$ with the mixed tensor $M^\mu{}_\nu=\bar g^{\mu\alpha}h_{\alpha\nu}$. Then $g^{-1}=(1+\epsilon M)^{-1}\bar g^{-1}$ and the Neumann series gives

$$
g^{\mu\nu}=\bar g^{\mu\nu}-\epsilon h^{\mu\nu}+\epsilon^2h^\mu{}_\rho h^{\rho\nu}-\epsilon^3h^\mu{}_\rho h^\rho{}_\sigma h^{\sigma\nu}+\cdots
$$

The implementation (`symbolic.metric_expansion.inverse_metric_series`) multiplies the truncated inverse by $g$ and checks that the product equals the identity through the declared order, for a generic diagonal background (four independent symbols) and generic symmetric $h$ (ten independent symbols), and also on a non-diagonal background [`exp.inverse_metric`, `exp.inverse_metric.nondiagonal`]. The explicit 2D example $\bar g={\rm diag}(-1,a^2)$ is an independent hand-checkable case [`exp.two_dim_example`].

## Determinant

For any $d$, $\det(1+\epsilon M)=\exp[{\rm tr}\ln(1+\epsilon M)]$, hence

$$
\frac{\sqrt{-g}}{\sqrt{-\bar g}}=\exp\Big[\tfrac12\epsilon\,{\rm tr}M-\tfrac14\epsilon^2{\rm tr}M^2+O(\epsilon^3)\Big]
=1+\frac\epsilon2h+\epsilon^2\Big(\frac18h^2-\frac14h_{\mu\nu}h^{\mu\nu}\Big)+O(\epsilon^3),
$$

using ${\rm tr}M=h$ and ${\rm tr}M^2=h^\mu{}_\nu h^\nu{}_\mu=h_{\mu\nu}h^{\mu\nu}$. The code does *not* use this formula to verify itself: it computes the exact determinant polynomial of $1+\epsilon M$ from the characteristic polynomial, takes the square root series, and compares [`exp.sqrt_det`]. Truncation errors scale as $\epsilon^{n+1}$ when the series is cut after order $n$, measured numerically (slopes 0.99, 2.01, 3.01 for the inverse and 0.97, 2.00, 3.01 for $\sqrt{-g}$, figure `geometry/expansion_truncation_errors.png`).

## Connection variation

The variation of the Levi-Civita connection is a tensor,

$$
\delta\Gamma^\lambda{}_{\mu\nu}=\tfrac12\bar g^{\lambda\rho}(\bar\nabla_\mu h_{\rho\nu}+\bar\nabla_\nu h_{\rho\mu}-\bar\nabla_\rho h_{\mu\nu}),
$$

derived from $\Gamma[\bar g+\epsilon h]$ by noting that the partial-derivative terms combine into covariant derivatives of $h$ because $\Gamma-\bar\Gamma$ is a tensor. It is verified against the direct derivative $\partial_\epsilon\Gamma$ at $\epsilon=0$ of the full connection of $\bar g+\epsilon h$ with a truncated inverse, for a generic polynomial $h$ [`lin.delta_gamma`].

## Second order

The second-order connection is $\delta^{(2)}\Gamma^\lambda{}_{\mu\nu}=-h^{\lambda\rho}(\bar\nabla_\mu h_{\rho\nu}+\bar\nabla_\nu h_{\rho\mu}-\bar\nabla_\rho h_{\mu\nu})/2$ (coefficient of $\epsilon^2$ in the exact connection, from the $-\epsilon h^{\lambda\rho}$ term of the inverse metric). First-order quantities are linear in $h$ and are (up to total derivatives) gauge-covariant objects built from $\bar\nabla h$. Second-order quantities contain products such as $h\,\bar\nabla\bar\nabla h$ and $(\bar\nabla h)^2$, and are what the quadratic action integrates (`docs/quadratic_action.md`). The tests exercise second order through the exact expansion of $\sqrt{-g}(R-2\Lambda)$ for the tensor ansatz rather than through a generic second-order Ricci tensor, so no statement is made here about the general second-order curvature beyond what that calculation supports.

## Sensitivity of the checks (negative controls)

| mistake | record | detected |
|---|---|---|
| wrong sign of $\epsilon^2$ term in $g^{-1}$ | `exp.inverse_metric.neg_sign` | yes |
| transposed index contraction (non-diagonal background) | `exp.inverse_metric.neg_index` | yes |
| factor of two in $h_{\mu\nu}h^{\mu\nu}$ term of $\sqrt{-g}$ | `exp.sqrt_det.neg_factor` | yes |

A transposed contraction is not detectable on a diagonal background with a symmetric $h$ in every case, which is why the control uses a non-diagonal background.
