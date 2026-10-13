# Data

`generated/` holds machine-generated tables. They are reproducible with `python scripts/run_parameter_sweeps.py` and `python scripts/benchmark_symbolics.py`; do not edit by hand.

| File | Content |
|---|---|
| sweep_wavenumber.csv | Spectrum versus comoving wavenumber $k$ |
| sweep_late_time.csv | Late-time freeze-out versus $x$ |
| sweep_H_Mpl.csv | Scaling with $H/M_{\rm pl}$, numeric vs predicted |
| sweep_bogoliubov.csv | Bogoliubov-rotated initial states, numeric vs predicted |
| sweep_x_start.csv | Effect of the start time $x_{\rm start}$ of leading-WKB initial data |
| sweep_tolerance.csv | Error versus integrator tolerance |
| symbolic_benchmark.csv | Timings of symbolic calculations |
| figure_catalog.json | Figure list with generating functions |

No external or observational data are used.
