# Gauge symmetry

## Diffeomorphisms and the sign convention

Let $\phi_\epsilon$ be the flow of a vector field $\xi^\mu$. The pulled-back metric is $\phi_\epsilon^*g=g+\epsilon\mathcal L_\xi g+O(\epsilon^2)$. A diffeomorphism changes the *description* of the same geometry; in coordinates $x'^\mu=x^\mu+\epsilon\xi^\mu$ the components as functions of the new coordinates are $g'_{\mu\nu}(x')=g_{\mu\nu}(x)-\epsilon(\mathcal L_\xi g)_{\mu\nu}$ (passive form). With $g=\bar g+\epsilon h$ and an unchanged background,

$$
h_{\mu\nu}\ \to\ h_{\mu\nu}-(\mathcal L_\xi\bar g)_{\mu\nu}=h_{\mu\nu}-\bar\nabla_\mu\xi_\nu-\bar\nabla_\nu\xi_\mu .
$$

The active convention corresponds to $\xi\to-\xi$; nothing physical depends on that choice, but every formula in the repository (gauge laws, invariants, residual freedom) uses the single choice above. The identity $\mathcal L_\xi\bar g=\bar\nabla\xi+\bar\nabla\xi$ (torsion-free, metric-compatible connection) is verified componentwise with generic polynomial $\xi$ [`gauge.lie_equals_covariant`].

## The field equation and gauge invariance

Because $\bar E_{\mu\nu}=0$ the linearized operator is gauge invariant:

$$
\delta E_{\mu\nu}[h-\mathcal L_\xi\bar g]=\delta E_{\mu\nu}[h]-\mathcal L_\xi\bar E_{\mu\nu}=\delta E_{\mu\nu}[h],
$$

verified for a generic $h$ [`gauge.field_equation_invariant`] and for pure gauge $h$ alone [`gauge.pure_gauge_solves_linearized`]. On a background that is not a solution the second term would not vanish and $\delta E$ would not be invariant; this is why the vacuum de Sitter background is used.

## Scalar-vector decomposition of the gauge parameter

Take conformal time and write the metric perturbation as

$$
ds^2=a^2\Big[-(1+2A)d\eta^2+2(\partial_iB+S_i)d\eta dx^i+\big((1-2\psi)\delta_{ij}+2\partial_i\partial_jE+\partial_iF_j+\partial_jF_i+\gamma_{ij}\big)dx^idx^j\Big],
$$

with $\partial_iS_i=\partial_iF_i=0$, $\partial_i\gamma_{ij}=\gamma_{ii}=0$. Decompose $\xi^\mu=(T,\partial^iL+L^i)$ with $\partial_iL^i=0$. Matching the components of $-\mathcal L_\xi\bar g$ (with $\bar g=a^2\eta_{\mu\nu}$, $\mathcal H=a'/a$) gives

$$
A\to A-(\mathcal HT+T'),\quad B\to B+T-L',\quad \psi\to\psi+\mathcal HT,\quad E\to E-L,
$$
$$
S_i\to S_i-L_i',\quad F_i\to F_i-L_i,\quad \gamma_{ij}\to\gamma_{ij}.
$$

These laws are obtained by symbolic component matching, not quoted [`gauge.svt_laws`]. Consequences:

* Only the *tensor* part $\gamma_{ij}$ is gauge invariant by itself.
* Combinations such as $\sigma=E'-B$ shift by $-T$, so $\Phi_B=A-\sigma'-\mathcal H\sigma$ and $\Psi_B=\psi+\mathcal H\sigma$ are invariant, and so is $V_i=S_i-F_i'$ (`docs/gauge_invariant_variables.md`).

## Gauge fixing and residual gauge freedom

| gauge | conditions | what remains |
|---|---|---|
| Newtonian | $B=E=0$ (and $S_i=F_i=0$) | none for $k\neq0$: $E\to E-L$ forces $L=0$, then $B\to B+T$ forces $T=0$ [`gauge.newtonian_no_residual`] |
| synchronous | $A=B=0$ | $T=C/a$ (so $aT$ constant) and $L'=T$: a function of space only [`gauge.synchronous_residual`] |
| spatially flat | $\psi=E=0$ | none: $\psi\to\psi+\mathcal HT$ with $\mathcal H=-1/\eta\neq0$ forces $T=0$, and $L=0$ |
| transverse-traceless (tensor) | $\partial_i\gamma_{ij}=\gamma_{ii}=0$ | none for the tensor sector itself (it is already gauge invariant) |

The $k=0$ sector is excluded from these statements: with $\partial_i\to0$ the functions $B,E,S_i,F_i$ do not exist as defined and the homogeneous parts must be treated separately (large gauge transformations; see `docs/limitations.md`).

## Constraints, gauge conditions and propagating degrees of freedom

* A **gauge condition** is imposed by hand to remove redundancy (4 conditions).
* A **constraint** is an equation without second time derivatives of some variable; the vacuum equations contain 4 of them ($\delta E_{00}$ and $\delta E_{0i}$), which determine non-dynamical metric components in terms of others.
* **Propagating degrees of freedom** are what remains: $10-4-4=2$, the two tensor polarizations (`docs/constraint_equations.md`).

## Numerical gauge experiment

`perturbations.gauge_invariants.numeric_gauge_experiment` builds one smooth scalar configuration and moves it through two non-trivial gauge shifts (analytic $T(\eta)$, $L(\eta)$, derivatives known exactly). The changes in $A,\psi,B,E$ are $O(10^{-1})$; the changes in $\Phi_B,\Psi_B$ are $\sim10^{-15}$ [`num.gauge_invariants`, `num.gauge_dependents`]. See figure `gauge/gauge_comparison.png`. The configuration is *kinematic* (not a vacuum solution): gauge invariance of $\Phi_B,\Psi_B$ does not rely on the field equations.

## Warning

A gauge-invariant expression is not automatically a local observable. $\Phi_B$ is invariant under *linearized* diffeomorphisms of the background, which is a statement about the perturbative description. See `docs/limitations.md`.
