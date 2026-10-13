# Verification Protocol

## Principle

A check that cannot fail proves nothing. Every verification is a `Record` (`validation/records.py`) with a name, a tier, a residual and a tolerance. Where practical a negative control is attached: the same check applied to a deliberately wrong input, which is expected to fail (`expect_fail`). A negative control that passes would flag a vacuous check.

## Tiers

1. Symbolic exact: residuals must simplify to zero (curvature of de Sitter, linearized Einstein tensor by two routes, Bianchi identities, gauge laws, quadratic action by two routes).
2. Analytic benchmarks: closed-form solutions (Bunch-Davies mode, power spectrum, Wronskian, Bogoliubov prediction).
3. Numerical: solver agreement across formulations, convergence orders, mpmath cross-check, SVT round-trips.

## Current status (python 3.13, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0)

- 128 of 128 records pass
- 9 negative controls behave as expected (they fail when they should)
- Runtime about 130 s, of which 118 s is symbolic

The machine-readable output is `reports/validation_report.json`; a readable version is `reports/validation_report.md`.

## Bugs caught by the protocol itself

During development the protocol exposed that sympy's `euler_equations` returned an empty list for a one-variable first-order check, making it vacuous. The project uses its own Euler-Lagrange operator, and a wrong-$\Lambda$ control now produces a non-zero residual.

## Running

    make validate        # full suite, writes reports/
    make test-fast       # pytest without slow marks
    make test            # everything
