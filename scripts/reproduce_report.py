#!/usr/bin/env python3
"""Reproduce everything: validation -> sweeps -> figures -> reports/generated_results.md.

Usage: python scripts/reproduce_report.py [--fast]
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(script, *args):
    print(">>", script, *args, flush=True)
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *args])
    return r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true")
    a = ap.parse_args()
    rc = run("run_validation_suite.py", *(["--fast"] if a.fast else []))
    run("run_parameter_sweeps.py")
    run("generate_all_figures.py")
    rep = json.loads((ROOT / "reports" / "validation_report.json").read_text())
    s = rep["summary"]
    by = {}
    for r in rep["records"]:
        by.setdefault(r["category"], []).append(r)
    lines = ["# Generated results", "", "Produced by `scripts/reproduce_report.py`; do not edit by hand.", "",
             f"Verification records: {s['total']} total, {s['passed']} passed, {s['failed']} failed "
             f"({s['negative_controls']} negative controls).", "",
             "| category | records | passed |", "|---|---|---|"]
    for c, rs in by.items():
        lines.append(f"| {c} | {len(rs)} | {sum(r['passed'] for r in rs)} |")
    key = [r for r in rep["records"] if r["id"] in (
        "num.exact_vs_conformal", "num.exact_vs_cosmic", "num.exact_vs_efold", "num.exact_vs_canonical",
        "num.bogoliubov_prediction", "num.spectrum_vs_exact", "num.rk4_order", "num.wronskian_conservation",
        "num.mpmath_crosscheck", "num.gauge_invariants")]
    lines += ["", "## Selected numerical residuals", "", "| id | residual | tolerance |", "|---|---|---|"]
    for r in key:
        lines.append(f"| `{r['id']}` | {r['residual']:.3e} | {r['tolerance']:.1e} |")
    (ROOT / "reports" / "generated_results.md").write_text("\n".join(lines) + "\n")
    print("wrote reports/generated_results.md")
    sys.exit(rc)


if __name__ == "__main__":
    main()
