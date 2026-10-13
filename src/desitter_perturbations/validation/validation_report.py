"""Run every check and write machine-readable (JSON) and human-readable (Markdown) reports."""

from __future__ import annotations

import json
import os
import platform
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import scipy
import sympy

from . import symbolic_checks, analytic_benchmarks, numerical_checks


def run_all(fast: bool = False, symbolic: bool = True):
    t0 = time.time()
    recs = []
    timings = {}
    if symbolic:
        t = time.time()
        recs += symbolic_checks.run_all()
        timings["symbolic"] = time.time() - t
    t = time.time()
    recs += analytic_benchmarks.run()
    timings["analytic"] = time.time() - t
    t = time.time()
    recs += numerical_checks.run(fast=fast)
    timings["numerical"] = time.time() - t
    meta = dict(
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sympy=sympy.__version__,
        platform=platform.platform(), seconds=round(time.time() - t0, 1), timings={k: round(v, 1) for k, v in timings.items()},
    )
    return recs, meta


def _cell(text):
    """Render a plain-text description as an inline code span that is safe inside a Markdown table."""
    return "`" + text.replace("`", "'").replace("|", "\\|") + "`"


def write_reports(recs, meta, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    passed = sum(r.passed for r in recs)
    summary = dict(total=len(recs), passed=passed, failed=len(recs) - passed,
                   negative_controls=sum(r.kind == "negative_control" for r in recs))
    payload = dict(summary=summary, environment=meta, records=[r.as_dict() for r in recs])
    (out / "validation_report.json").write_text(json.dumps(payload, indent=2))
    by_cat = defaultdict(list)
    for r in recs:
        by_cat[r.category].append(r)
    lines = ["# Validation report", "",
             f"Checks: {summary['total']}; passed: {summary['passed']}; failed: {summary['failed']}; "
             f"of which negative controls: {summary['negative_controls']}.", "",
             "Environment: " + ", ".join(f"{k} {v}" for k, v in meta.items() if k not in ("timings",)), "",
             "Symbolic residuals are 0 (identically zero after simplification) or 1. A *negative control* passes only when the deliberately wrong "
             "formula is detected, i.e. its residual exceeds the tolerance.", ""]
    for cat, rs in by_cat.items():
        lines += [f"## {cat}", "", "| id | kind | description | residual | tolerance | result |", "|---|---|---|---|---|---|"]
        for r in rs:
            lines.append(f"| `{r.id}` | {r.kind} | {_cell(r.description)} | {r.residual:.3e} | {r.tolerance:.1e} | {'PASS' if r.passed else 'FAIL'} |")
        lines.append("")
    (out / "validation_report.md").write_text("\n".join(lines))
    return summary
