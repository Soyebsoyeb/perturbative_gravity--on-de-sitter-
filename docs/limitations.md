# Limitations and Scope

## What this project does

- Derives and verifies linearized gravity on exact de Sitter, gauge structure, SVT decomposition, the tensor sector action, free canonical quantization of gravitons and the resulting spectrum.

## What it does not do

- It does not solve or claim progress on quantum gravity. The quantization is of free quadratic fluctuations on a fixed background.
- No interacting vertices, loops, backreaction, or non-linear gauge fixing.
- Scalar and vector sectors are analysed classically (constraints, gauge invariants), with no matter content. Pure de Sitter with a cosmological constant has no propagating scalar or vector gravity modes; the project does not treat inflaton fluctuations.
- Exact de Sitter only. Slow-roll corrections, spectral tilt and transfer functions are not included.
- The Bunch-Davies state is assumed as the physical vacuum; its choice is a physical input, not a derivation.
- Results depend on conventions (signature, polarization normalization). They are fixed in `conventions.md` and checked, but external comparisons must convert.

## Numerical caveats

- Residual diagnostics measure interpolation error as well as solver error.
- Precision is limited by tolerance settings; double-precision floors are near 1e-10 for the cross-formulation comparison.
- The periodic-box SVT decomposition has zero-mode ambiguities handled by an explicit rule.

## Exploratory content

The finite-$x_{\rm start}$ Bogoliubov interpretation of leading-WKB initial data is an exact consequence of the equations, but its presentation as a systematic is exploratory; no claim of novelty beyond a pedagogical clarification is made.

## Truncated specification

The master build prompt was cut off in section 11 (direct tensor evolution and canonical quantization). The remaining items were implemented from the visible requirements and the project principles.
