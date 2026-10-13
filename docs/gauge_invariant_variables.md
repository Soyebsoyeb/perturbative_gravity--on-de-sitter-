# Gauge-invariant variables

## Definitions

Using the transformation laws of `docs/gauge_symmetry.md` ($\sigma\equiv E'-B$, $\sigma\to\sigma-T$):

$$
\Phi_B=A-\frac1a(a\sigma)'=A-\sigma'-\mathcal H\sigma,\qquad
\Psi_B=\psi+\mathcal H\sigma,\qquad
V_i=S_i-F_i',\qquad \gamma_{ij}.
$$

Check: $\Phi_B\to A-(\mathcal HT+T')-(\sigma'-T')-\mathcal H(\sigma-T)=\Phi_B$; $\Psi_B\to\psi+\mathcal HT+\mathcal H(\sigma-T)=\Psi_B$; $V_i\to S_i-L_i'-F_i'+L_i'=V_i$. Verified symbolically (gauge laws) and numerically for two independent gauge shifts at the $10^{-15}$ level [`num.gauge_invariants`].

## What they do

* $\Phi_B$ and $\Psi_B$ are the lapse and curvature perturbations in Newtonian gauge ($B=E=0$): in that gauge $\Phi_B=A$, $\Psi_B=\psi$.
* All components of the linearized field equations depend on the metric perturbation only through these variables (scalar sector) or $V_i$ (vector sector) or $\gamma_{ij}$, verified exactly [`sector.*`].
* The vacuum constraints force $\Phi_B=\Psi_B=0$ and $V_i=0$ for $k\neq0$ (`docs/constraint_equations.md`).

## Different gauges, same physical perturbation

The comparison in `figures/gauge/gauge_comparison.png` shows $A$ and $\psi$ in three gauges differing by $O(1)$ and diverging near $\eta\to0$ where $\mathcal HT$ grows, while $\Phi_B$ and $\Psi_B$ coincide. The coordinate-dependent components therefore carry no meaning on their own.

## Limits of the concept

1. **Invariance is perturbative.** These quantities are invariant under *linearized* gauge transformations. Beyond linear order the transformation laws acquire terms quadratic in the fields and in $\xi$, and the invariants must be redefined order by order.
2. **Invariant does not mean observable.** In a diffeomorphism-invariant theory, local observables require a physical reference frame (relational observables built from matter or boundary data). $\Phi_B$ depends on the background slicing used in its definition (the choice of $\eta$, i.e. a preferred foliation) and is not a local scalar built from curvature.
3. **Vacuum de Sitter.** Here $\Phi_B$ vanishes anyway; its role is to show *how* the dynamics reduces to the tensor sector.
4. **Tensor sector.** $\gamma_{ij}$ is invariant under the linearized transformations only in the sense that its transformation law is trivial; at higher order its definition is gauge dependent.
