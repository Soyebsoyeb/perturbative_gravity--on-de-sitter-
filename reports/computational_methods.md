# Computational Methods

## Symbolic engine

Pure SymPy. Geometry is computed by a `Geometry` class parameterized by the Riemann sign, so convention tests can be run with the wrong sign. Perturbation expansions use a bookkeeping parameter $\epsilon$ and truncate by series. Heavy results are cached with `lru_cache`. Typical timings: SVT algebra 32 s, direct second-order action 3.8 s, route-two action 2.3 s; see `data/generated/symbolic_benchmark.csv`.

## Verification engine

`Record` objects with a kind (identity, benchmark or negative control), residual, tolerance. Symbolic residuals reduce to 0 or 1; numerical residuals are relative errors. Reports are written in JSON and Markdown by `validation_report.write_reports`.

## Numerical engine

Four ODE formulations in a common state layout; DOP853 at $\mathrm{rtol}=10^{-10}$ by default; logarithmic output in $x=k/(aH)$; fixed-step RK4 in e-folds for order measurement; mpmath odefun for high precision.

## Spectral SVT decomposition

Periodic FFT box, spectral and finite-difference divergences, explicit zero-mode rule $\psi_0=-T/6$, noise amplification measured for $E$ and $F_i$.

## Plotting

Matplotlib with a fixed categorical palette order, recessive grid, direct labels and at most four series per panel. Figures were inspected and corrected for clipped labels, legend overlap and artifacts.

## Reproducibility

Deterministic seeds, pinned numbers recorded in the validation report environment block, and `make all` rebuilding results end to end.
