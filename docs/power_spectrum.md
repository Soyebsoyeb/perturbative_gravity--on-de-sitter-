# Tensor Power Spectrum

Status: exact analytic result, reproduced numerically to about $10^{-11}$.

## Definition

With the polarization normalization $e_{ij}e_{ij}=2$ and two polarizations,

$$
\Delta_t^2(k)=\frac{k^3}{2\pi^2}\sum_{\lambda}2\,\lvert\gamma_k\rvert^2=\frac{2k^3}{\pi^2}\lvert\gamma_k\rvert^2 ,
$$

where $\gamma_k$ is the mode function of one polarization amplitude, so that $\langle\gamma_{ij}\gamma_{ij}\rangle=\int d\ln k\,\Delta_t^2$.

## Exact de Sitter result

With the Bunch-Davies mode,

$$
\Delta_t^2(k,\eta)=\frac{2H^2}{\pi^2M_{\rm pl}^2}\,\big(1+x^2\big),\qquad x=\frac{k}{aH}=-k\eta .
$$

At late times $x\to0$ this freezes to the scale-invariant value $2H^2/(\pi^2M_{\rm pl}^2)$. At early times ($x\gg1$) the expression grows as $x^2$ because the mode is still oscillating inside the horizon.

## Numerical agreement

The ODE-based spectrum matches the exact expression with maximum relative error $9.9\times10^{-12}$ over the output grid (`validation/analytic_benchmarks.py`).

## Dependence on the initial state

Leading-WKB initial data at finite $x_{\rm start}$ is a Bogoliubov state with $\lvert\beta\rvert=1/(2x_{\rm start}^2)$. The late-time power is multiplied by $\lvert\alpha-\beta\rvert^2$, and the numerical value agrees with this prediction (`data/generated/sweep_x_start.csv`).

## Scaling

$\Delta_t^2\propto(H/M_{\rm pl})^2$ (`sweep_H_Mpl.csv`) and is independent of $k$ at late times (`sweep_wavenumber.csv`).

## Caveat

This is the spectrum of a free tensor field on exact de Sitter. A realistic inflationary spectrum has a small tilt from the slow-roll evolution of the background, which is outside the present scope.
