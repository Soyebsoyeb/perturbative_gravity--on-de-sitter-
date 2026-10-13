# Numerical Methods

## Formulations

The tensor mode equation is integrated in four independent formulations, each with a real four-component state (real and imaginary parts of the value and its derivative):

1. conformal time: $\gamma''+2\mathcal H\gamma'+k^2\gamma=0$
2. cosmic time
3. e-fold number
4. canonical variable: $v''+\big(k^2-2/\eta^2\big)v=0$

The default integrator is `scipy.integrate.solve_ivp` with DOP853 and $\mathrm{rtol}=10^{-10}$.

## Results at $\mathrm{rtol}=10^{-10}$

| Formulation | max relative error vs exact |
|---|---|
| conformal | $7.6\times10^{-10}$ |
| cosmic | $5.5\times10^{-10}$ |
| e-fold | $5.5\times10^{-10}$ |
| canonical | $7.9\times10^{-10}$ |
| cross-formulation | about $3.1\times10^{-10}$ |

## Convergence

A fixed-step RK4 in e-folds has observed order $4.003$. The tolerance sweep (`sweep_tolerance.csv`) shows the error tracking the requested tolerance until it saturates at the floor set by the conversion between formulations.

## Independent cross-check

`mpmath.odefun` at high precision agrees with the double-precision solution to $2.6\times10^{-13}$ (slow test).

## Residual diagnostics

The mode is interpolated with a cubic spline and substituted into the ODE. The residual converges at order $4.00$ in output density, with a peak of $2.4\times10^{-6}$ at 4000 points. It diagnoses the interpolant, not the integrator.

## Output grid

Output is on a logarithmic grid in $x=k/(aH)$, which resolves horizon crossing; a linear grid in time would put almost all points deep inside the horizon.

## SVT numerics

The scalar-vector-tensor decomposition on a periodic FFT box uses spectral derivatives. Zero modes follow the rule $\psi_0=-T/6$. Small-$k$ amplification of the reconstructed potentials ($E\sim k^{-2}$, $F\sim k^{-1}$) is measured in `noise_amplification`.
