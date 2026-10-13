#!/usr/bin/env python3
"""Run the parameter sweeps and write CSV tables to data/generated/."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from desitter_perturbations.spectra import parameter_sweeps as ps  # noqa: E402

OUT = ROOT / "data" / "generated"


def main():
    tables = {
        "sweep_wavenumber": ps.sweep_wavenumber(n=30, k_min=0.1, k_max=200.0, eta_end=-1e-2, x_start=300.0),
        "sweep_late_time": ps.sweep_late_time(),
        "sweep_H_Mpl": ps.sweep_H_Mpl(),
        "sweep_bogoliubov": ps.sweep_bogoliubov(),
        "sweep_x_start": ps.sweep_x_start(),
        "sweep_tolerance": ps.sweep_tolerance(),
    }
    for name, t in tables.items():
        ps.write_csv(str(OUT / f"{name}.csv"), t)
        print("wrote", OUT / f"{name}.csv", f"({len(next(iter(t.values())))} rows)")


if __name__ == "__main__":
    main()
