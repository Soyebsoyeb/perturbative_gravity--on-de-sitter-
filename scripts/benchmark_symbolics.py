#!/usr/bin/env python3
"""Time the principal symbolic computations (written to data/generated/symbolic_benchmark.csv)."""
import csv
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import sympy as sp  # noqa: E402
from desitter_perturbations import background as bg  # noqa: E402
from desitter_perturbations.symbolic import curvature as cu, einstein_tensor as et, quadratic_action as qa  # noqa: E402
from desitter_perturbations.symbolic import tensor_decomposition as td  # noqa: E402
from desitter_perturbations.perturbations import scalar_sector  # noqa: E402


def timed(label, fn, rows):
    t = time.perf_counter()
    fn()
    dt = time.perf_counter() - t
    rows.append((label, round(dt, 3)))
    print(f"{label:55s} {dt:8.2f} s")


def main():
    b = bg.symbolic_background()
    rows = []
    g = cu.Geometry(b["coords_conformal"], b["g_conformal"])
    timed("de Sitter Riemann + invariants (conformal)", lambda: (g.riemann_down(), g.kretschmann()), rows)
    h = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f"h{min(i,j)}{max(i,j)}") * b["eta"])
    timed("linearized Einstein, covariant (generic symbols)", lambda: et.linearized_einstein_covariant(g, h, 3 * b["H"] ** 2), rows)
    timed("quadratic action, route 1 (direct expansion)", lambda: qa.route1_direct_expansion.__wrapped__(), rows)
    timed("quadratic action, route 2 (linearized equations)", lambda: qa.route2_from_linearized_equations.__wrapped__(), rows)
    timed("SVT algebra identities (generic h, k)", lambda: td.verify_svt_identities.__wrapped__(), rows)
    timed("scalar sector in general gauge", lambda: scalar_sector.general_gauge_dE.__wrapped__(), rows)
    out = ROOT / "data" / "generated" / "symbolic_benchmark.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["computation", "seconds"])
        w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
