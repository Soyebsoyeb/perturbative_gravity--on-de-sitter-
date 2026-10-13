# Tensor-mode dynamics in de Sitter space

## Equations

From the quadratic action (`docs/quadratic_action.md`) and independently from $\delta E_{ij}$, each polarization obeys

$$
\gamma_k''+2\mathcal H\gamma_k'+k^2\gamma_k=0,\qquad\mathcal H=-1/\eta,
$$

equivalently in cosmic time $\ddot\gamma_k+3H\dot\gamma_k+(k^2/a^2)\gamma_k=0$ (using $\gamma'=a\dot\gamma$, $\gamma''=a^2\ddot\gamma+a^2H\dot\gamma$ so that $\gamma''+2\mathcal H\gamma'=a^2(\ddot\gamma+3H\dot\gamma)$; the identity is tested in `test_cosmic_and_conformal_equations_equivalent`), and in e-folds $N=\ln a$: $\gamma_{NN}+3\gamma_N+(k/(aH))^2\gamma=0$.

## Exact solutions

Try $\gamma=(1+ik\eta)e^{-ik\eta}$: $\gamma'=k^2\eta\,e^{-ik\eta}$, $\gamma''=(k^2-ik^3\eta)e^{-ik\eta}$, so

$$
\gamma''-\frac2\eta\gamma'+k^2\gamma=e^{-ik\eta}\big[k^2-ik^3\eta-2k^2+k^2+ik^3\eta\big]=0 .
$$

The complex conjugate is the second solution. The normalised Bunch-Davies mode (`docs/quantum_state.md`) is

$$
\gamma_k(\eta)=\frac{iH}{M_{\rm pl}k^{3/2}}(1+ik\eta)e^{-ik\eta},\qquad |\gamma_k|^2=\frac{H^2}{M_{\rm pl}^2k^3}(1+x^2),\ x=-k\eta .
$$

The residual of this solution in the ODE is $\lesssim7\times10^{-16}$ over $k\in[0.05,50]$, $|\eta|\in[10^{-3},10^3]$ [`bench.mode_solves_ode`].

## Asymptotic regimes

| regime | behaviour | benchmark |
|---|---|---|
| subhorizon, $x\gg1$ | $\gamma\to e^{-ik\eta}/(aM_{\rm pl}\sqrt k)$, amplitude $\propto1/a$; relative error of $\lvert\gamma\rvert$ is $1/(2x^2)$ | `bench.subhorizon.*` |
| horizon crossing, $x=1$ | $\lvert\gamma\rvert=\sqrt2\,H/(M_{\rm pl}k^{3/2})$ exactly | `test_horizon_crossing_amplitude_is_sqrt2` |
| superhorizon, $x\ll1$ | $\gamma\to iH/(M_{\rm pl}k^{3/2})$ (frozen); relative correction $x^2/2$ | `bench.superhorizon.*` |

The subhorizon form applies when the physical wavelength is much smaller than the Hubble radius; the freeze-out follows because the friction term $2\mathcal H\gamma'$ dominates over $k^2\gamma$ once $k<\mathcal H$.

## Three viewpoints on a mode

* **Conformal time**: $\eta$ is the natural variable because the equation is that of a time-dependent oscillator; the comoving horizon $-\eta$ shrinks, so modes *exit* the horizon.
* **Cosmic time**: $\gamma$ is damped by Hubble friction $3H\dot\gamma$ and its gradient term redshifts as $k^2e^{-2Ht}$.
* **Physical wavelength**: $\lambda_{\rm phys}=2\pi a/k$ grows exponentially while the Hubble radius $1/H$ is constant; crossing occurs at $\lambda_{\rm phys}=2\pi/H$.

## Canonical variable

With $z_T=aM_{\rm pl}/\sqrt2$ (derived in `docs/canonical_normalization.md`), $v=z_T\gamma$ obeys

$$
v''+\Big(k^2-\frac{z_T''}{z_T}\Big)v=0,\qquad \frac{z_T''}{z_T}=\frac{a''}{a}=\frac2{\eta^2}.
$$

This is the equation of a massless, minimally coupled *scalar* field rescaled by $a$. The resemblance is a mathematical analogy: tensor perturbations and scalar fields have different origins, normalisations ($z_T$ versus $a$) and interpretations. The solutions coincide because $z_T\propto a$ in exact de Sitter space.

## Numerical solutions

Four formulations (conformal, cosmic, e-fold, canonical) agree with the exact solution to $5$--$8\times10^{-10}$ at $\mathrm{rtol}=10^{-10}$ and with each other to $3\times10^{-10}$ [`num.exact_vs_*`, `num.cross_*`]; see `docs/numerical_methods.md`.
