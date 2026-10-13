# Linearized curvature and the linearized Einstein equation

## Covariant formulas

With $\delta\Gamma$ from `docs/metric_variation.md`, the Palatini identity gives

$$
\delta R_{\mu\nu}=\bar\nabla_\lambda\delta\Gamma^\lambda{}_{\mu\nu}-\bar\nabla_\nu\delta\Gamma^\lambda{}_{\lambda\mu},
\qquad
\delta R=-h^{\mu\nu}\bar R_{\mu\nu}+\bar g^{\mu\nu}\delta R_{\mu\nu},
$$
$$
\delta G_{\mu\nu}=\delta R_{\mu\nu}-\tfrac12h_{\mu\nu}\bar R-\tfrac12\bar g_{\mu\nu}\delta R,\qquad
\delta E_{\mu\nu}=\delta G_{\mu\nu}+\Lambda h_{\mu\nu}.
$$

Expanding $E_{\mu\nu}=\bar E_{\mu\nu}+\epsilon\,\delta E_{\mu\nu}+O(\epsilon^2)$ with $E=G+\Lambda g$ and $g=\bar g+\epsilon h$ yields $\delta E=\delta G+\Lambda h$ directly, because $\delta g_{\mu\nu}=h_{\mu\nu}$ for the covariant metric.

## Simplification on de Sitter space

On the background solution $\bar R_{\mu\nu}=3H^2\bar g_{\mu\nu}$, $\bar R=12H^2$, $\Lambda=3H^2$, and the explicit curvature terms combine:

$$
\delta E_{\mu\nu}=\delta R_{\mu\nu}-3H^2h_{\mu\nu}+\tfrac32H^2\bar g_{\mu\nu}h-\tfrac12\bar g_{\mu\nu}\big(\bar\nabla_\alpha\bar\nabla_\beta h^{\alpha\beta}-\bar\Box h\big),
$$

where $\delta R_{\mu\nu}=\tfrac12(\bar\nabla^\lambda\bar\nabla_\mu h_{\nu\lambda}+\bar\nabla^\lambda\bar\nabla_\nu h_{\mu\lambda}-\bar\Box h_{\mu\nu}-\bar\nabla_\mu\bar\nabla_\nu h)$. The cosmological constant therefore enters twice: through the background curvature in $\delta R$ and $\delta G$ (via $\bar R_{\mu\nu}$, $\bar R$, and via the commutator of covariant derivatives in $\delta R_{\mu\nu}$) and explicitly through $\Lambda h_{\mu\nu}$. Terms that survive *because* the curvature is constant: with $\bar R_{\mu\nu\rho\sigma}=H^2(\bar g\bar g-\bar g\bar g)$ every commutator $[\bar\nabla,\bar\nabla]h$ becomes algebraic in $h$ with coefficients $H^2$, and no derivative of the background curvature ever appears.

For a transverse ($\bar\nabla^\mu h_{\mu\nu}=0$), traceless perturbation this reduces to the massless spin-2 equation

$$
\delta E_{\mu\nu}=-\tfrac12\big(\bar\Box-2H^2\big)h_{\mu\nu}.
$$

This statement is verified exactly for the plus/cross ansatz: the ansatz is transverse and traceless in the covariant sense, and the engine's $\delta E_{\mu\nu}$ equals $-\tfrac12(\bar\Box h_{\mu\nu}-2H^2h_{\mu\nu})$ with $\bar\Box$ computed from two componentwise covariant derivatives [`sector.tensor.covariant`].

## Two independent derivations

* **Covariant** (`einstein_tensor.linearized_einstein_covariant`): the formulas above evaluated componentwise.
* **Component level** (`einstein_tensor.linearized_einstein_bruteforce`): build $g=\bar g+\epsilon h$, compute a truncated inverse, the full Christoffel symbols and Ricci tensor, truncate at $O(\epsilon)$ after every product, and extract the coefficient of $\epsilon$ of $G+\Lambda g$. No covariant derivative of $h$ occurs.

They agree identically (exact rational-function cancellation) for a generic polynomial $h$ with all ten components non-zero, in conformal coordinates and again in cosmic-time coordinates [`lin.cov_vs_component`, `lin.cov_vs_component.cosmic`]. The brute-force pipeline also returns $\bar E_{\mu\nu}=0$ at $O(\epsilon^0)$, i.e. the background field equations [`lin.background_zero`].

## Bianchi identity

The exact identity $\nabla_\mu G^\mu{}_\nu=0$ linearizes (using $\bar\nabla\bar G^\mu{}_\nu=0$) to

$$
\bar\nabla_\mu\delta G^\mu{}_\nu+\delta\Gamma^\mu{}_{\mu\lambda}\bar G^\lambda{}_\nu-\delta\Gamma^\lambda{}_{\mu\nu}\bar G^\mu{}_\lambda=0,\qquad
\delta G^\mu{}_\nu=\bar g^{\mu\rho}\delta G_{\rho\nu}-h^{\mu\rho}\bar G_{\rho\nu}.
$$

In covariant-index form the same statement is

$$
\bar g^{\mu\rho}\big[\bar\nabla_\rho\delta G_{\mu\nu}-\delta\Gamma^\lambda{}_{\rho\mu}\bar G_{\lambda\nu}-\delta\Gamma^\lambda{}_{\rho\nu}\bar G_{\mu\lambda}\big]=0 .
$$

The naive divergence $\bar\nabla^\mu\delta G_{\mu\nu}$ is **not** zero off the field equations: the remainder is proportional to $\bar G_{\mu\nu}=-\Lambda\bar g_{\mu\nu}$ and is displayed in notebook 03. For $E=G+\Lambda g$ the background value $\bar E=0$ removes those terms and $\bar\nabla^\mu\delta E_{\mu\nu}=0$ exactly. All three statements, and the non-vanishing of the naive divergence (a negative control), are verified [`bianchi.covariant_form`, `bianchi.mixed_form`, `bianchi.field_equation_divergence`, `bianchi.naive_divergence_nonzero`].

## Pure-gauge perturbations

Since $\bar E=0$ and $\mathcal L_\xi\bar E=0$, a pure-gauge $h=-\mathcal L_\xi\bar g$ solves $\delta E=0$, and the linearized operator is gauge invariant: $\delta E[h+h_{\rm pg}]=\delta E[h]$ [`gauge.pure_gauge_solves_linearized`, `gauge.field_equation_invariant`].
