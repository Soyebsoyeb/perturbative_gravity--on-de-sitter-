# Second-order Einstein-Hilbert action

## Statement

For the transverse-traceless tensor ansatz

$$
ds^2=a^2\big[-d\eta^2+(\delta_{ij}+\epsilon\gamma_{ij})dx^idx^j\big],\quad \gamma_{xx}=-\gamma_{yy}=f_+(\eta,z),\ \gamma_{xy}=f_\times(\eta,z),
$$

the second-order Lagrangian density (coefficient of $1/16\pi G$) is, up to total derivatives,

$$
\mathcal L_2=\frac{a^2}{4}\big(\gamma_{ij}'\gamma_{ij}'-\partial_k\gamma_{ij}\partial_k\gamma_{ij}\big)
\ \Rightarrow\ S_2=\frac{M_{\rm pl}^2}{8}\int d\eta\,d^3x\,a^2\big[\gamma_{ij}'\gamma_{ij}'-(\partial_k\gamma_{ij})^2\big].
$$

Here $\gamma_{ij}\gamma_{ij}=2(f_+^2+f_\times^2)$ and $1/(16\pi G)=M_{\rm pl}^2/2$. This is **derived and verified** here; the quadratic action is never copied.

## Expansion bookkeeping

$S=\frac1{16\pi G}\int d^4x\,\mathcal L$, $\mathcal L=\sqrt{-g}(R-2\Lambda)$. Exact expansion:

| order | result | note |
|---|---|---|
| $\epsilon^0$ | $\mathcal L_0=6/(H^2\eta^4)$ | background action density [`action.background_lagrangian`] |
| $\epsilon^1$ | $0$ for the TT ansatz; for a generic diagonal $h(\eta)$ the term is a total derivative | [`action.first_order_tensor`, `action.first_order_total_derivative`] |
| $\epsilon^2$ | $\mathcal L_2$ | below |

The first-order check is not vacuous: with $\Lambda=2H^2$ instead of $3H^2$ the first-order term is **not** a total derivative [`action.first_order_wrong_lambda`], which is the statement that the background does not solve the field equations. The first-order term equals $-\sqrt{-\bar g}\bar E^{\mu\nu}h_{\mu\nu}$ plus a total derivative, and vanishes on the solution.

## Route 1: direct expansion

Compute the exact metric inverse (block structure), $\sqrt{-g}=a^4\sqrt{1-\epsilon^2(f_+^2+f_\times^2)}$, the exact Ricci scalar, and the Taylor coefficient of $\epsilon^2$. The result contains second derivatives ($f\,f''$) and terms such as $\eta f f'$ and $f^2$; integrating by parts turns it into $\mathcal L_2$. The difference $\mathcal L_2^{\rm direct}-\mathcal L_2$ is a total derivative iff its Euler-Lagrange expressions vanish identically, which is what is checked [`action.route1_vs_target`].

The integrations by parts implicit in this equality are:

1. $-4\eta^2 f(f''-f_{zz})\to+4\eta^2(f'^2-f_z^2)+\partial_\eta(\cdot)+\partial_z(\cdot)$, with boundary terms $\partial_\eta(-4\eta^2ff')$ and $\partial_z(4\eta^2ff_z)$;
2. $12\eta ff'=6\eta\,\partial_\eta(f^2)\to-6f^2+\partial_\eta(6\eta f^2)$, boundary term $6\eta f^2$.

(The coefficients are those of the verified expression; all boundary terms are total derivatives and are discarded.) The spatial boundary terms vanish for periodic or compactly supported $\gamma$; the time-boundary terms are the standard Gibbons-Hawking-York-type terms that fix the variational problem, and they do not affect the equations of motion.

## Route 2: from the linearized equations

Because $\delta S=-\frac1{16\pi G}\int\sqrt{-g}\,E^{\mu\nu}\delta g_{\mu\nu}$ and $\bar E=0$, the second variation is

$$
\mathcal L_2=-\tfrac12\sqrt{-\bar g}\,h^{\mu\nu}\delta E_{\mu\nu}[h],
$$

with $\delta E$ from the covariant linearization of `docs/linearized_curvature.md`. This is an independent calculation (no second-order curvature is used) and also equals $\mathcal L_2$ up to a total derivative [`action.route2_vs_target`]. The two routes agree with each other [`action.route1_vs_route2`]. A doubled normalisation is rejected by the same test [`action.normalisation_control`].

## Equation of motion

The Euler-Lagrange equation of $\mathcal L_2$ for either polarization is $-a^2(\gamma''+2\mathcal H\gamma'-\nabla^2\gamma)=0$ [`action.tensor_eom`], identical to the equation from $\delta E_{ij}=\tfrac12(\gamma''+2\mathcal H\gamma'-\nabla^2\gamma)$, so the action and the equation of motion use the same $\gamma_{ij}$, the same $M_{\rm pl}$ and the same polarization normalisation.

## Field redefinitions

Two parametrisations are common: $g_{ij}=a^2(\delta_{ij}+\gamma_{ij})$ (used here) and $g_{ij}=a^2(e^{\gamma})_{ij}$. They differ by $\gamma_{ij}\to\gamma_{ij}+\tfrac12\gamma_{ik}\gamma_{kj}$. The $O(\gamma^2)$ term shifts the second-order action by $\tfrac12(\gamma\gamma)_{ij}\,\delta S/\delta g_{ij}\big|_{\rm bg}$, which vanishes because the background satisfies the field equation. The quadratic actions therefore agree up to boundary terms. At cubic order the parametrisations differ and the results are not interchangeable (not used here).

## What is not done

The general (non-TT) second-order action, with scalar and vector parts, was not derived. Gauge-fixing the full quadratic action and solving the constraints is outside this repository's implemented scope; the constraint analysis of `docs/constraint_equations.md` is done at the level of the linearized equations.
