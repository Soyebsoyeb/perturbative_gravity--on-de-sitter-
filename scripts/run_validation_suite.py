#!/usr/bin/env python3
"""Run every symbolic, analytic and numerical check; write reports/validation_report.{json,md}.

Usage: python scripts/run_validation_suite.py [--fast] [--no-symbolic]
Exit status is non-zero if any check fails.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from desitter_perturbations.validation import validation_report as vr  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="skip the arbitrary-precision cross-check")
    ap.add_argument("--no-symbolic", action="store_true", help="skip the (slow) symbolic checks")
    ap.add_argument("--out", default=str(ROOT / "reports"))
    a = ap.parse_args()
    recs, meta = vr.run_all(fast=a.fast, symbolic=not a.no_symbolic)
    summary = vr.write_reports(recs, meta, a.out)
    print(f"{summary['passed']}/{summary['total']} checks passed "
          f"({summary['negative_controls']} negative controls); reports in {a.out}")
    for r in recs:
        if not r.passed:
            print("FAILED:", r.id, r.residual, r.tolerance)
    sys.exit(0 if summary["failed"] == 0 else 1)


if __name__ == "__main__":
    main()
