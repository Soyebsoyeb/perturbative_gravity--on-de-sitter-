import json
import os
import numpy as np

from desitter_perturbations.perturbations import svt_numeric as sv, gauge_invariants as gi
from desitter_perturbations.spectra import parameter_sweeps as ps
from desitter_perturbations.validation import analytic_benchmarks, numerical_checks, validation_report
from desitter_perturbations.validation.records import record


def test_random_fields_are_seeded():
    assert np.array_equal(sv.random_symmetric_field(8, seed=3), sv.random_symmetric_field(8, seed=3))
    assert not np.array_equal(sv.random_symmetric_field(8, seed=3), sv.random_symmetric_field(8, seed=4))


def test_numerics_are_deterministic():
    a = ps.sweep_wavenumber(n=5, k_min=0.3, k_max=20.0, x_start=100.0)
    b = ps.sweep_wavenumber(n=5, k_min=0.3, k_max=20.0, x_start=100.0)
    assert np.array_equal(a["numeric"], b["numeric"])
    e1 = gi.numeric_gauge_experiment()[1]["shifted_1"]["Phi"]
    e2 = gi.numeric_gauge_experiment()[1]["shifted_1"]["Phi"]
    assert np.array_equal(e1, e2)


def test_figure_generation_is_deterministic(tmp_path, monkeypatch):
    from desitter_perturbations.visualization import spectrum_plots
    monkeypatch.setenv("DESITTER_FIGDIR", str(tmp_path / "a"))
    p1 = spectrum_plots.plot_scaling_H_M()
    monkeypatch.setenv("DESITTER_FIGDIR", str(tmp_path / "b"))
    p2 = spectrum_plots.plot_scaling_H_M()
    assert open(p1, "rb").read() == open(p2, "rb").read()


def test_report_roundtrip(tmp_path):
    recs = analytic_benchmarks.run()
    summary = validation_report.write_reports(recs, {"note": "test"}, tmp_path)
    data = json.loads((tmp_path / "validation_report.json").read_text())
    assert data["summary"]["total"] == len(recs) == summary["total"]
    assert data["summary"]["failed"] == 0
    assert all({"id", "residual", "tolerance", "passed", "assumptions"} <= set(r) for r in data["records"])
    assert (tmp_path / "validation_report.md").read_text().startswith("# Validation report")


def test_negative_control_semantics():
    ok = record("x", "c", "d", residual=1.0, tolerance=0.0, expect_fail=True)
    bad = record("y", "c", "d", residual=0.0, tolerance=0.0, expect_fail=True)
    assert ok.passed and not bad.passed


def test_fast_numerical_suite_passes():
    recs = numerical_checks.run(fast=True)
    failed = [r.id for r in recs if not r.passed]
    assert not failed, failed
