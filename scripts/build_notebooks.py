#!/usr/bin/env python3
"""Build and execute the 14 study notebooks (outputs are embedded).

Usage: python scripts/build_notebooks.py [--no-execute] [--only 07 08]
"""
import argparse
import sys
from pathlib import Path

import nbformat as nbf
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks"

SETUP = '''import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "src")))
import numpy as np, sympy as sp
from IPython.display import Image, display, Markdown
sp.init_printing(use_unicode=True)
def show(path):
    p = path[0] if isinstance(path, tuple) else path
    display(Image(filename=p))
'''

def md(t): return ("md", t.strip())
def code(t): return ("code", t.strip("\n"))

NOTEBOOKS = {}

NOTEBOOKS["01_de_sitter_geometry"] = ("Geometry of flat-sliced de Sitter space", [
md("""
# 01. Geometry of flat-sliced de Sitter space

Goal: build the background in both coordinate systems, derive the relation between cosmic and conformal time,
and verify -- by direct computation with the conventions of `docs/conventions.md` -- the constant-curvature identity
$R_{\\mu\\nu\\rho\\sigma}=H^2(g_{\\mu\\rho}g_{\\nu\\sigma}-g_{\\mu\\sigma}g_{\\nu\\rho})$, the Einstein tensor and the curvature invariants.
Every sign is checked against the implementation, including a negative control with a deliberately flipped Riemann sign.
"""),
code("""
from desitter_perturbations import background as bg
from desitter_perturbations.symbolic.curvature import Geometry
from desitter_perturbations.symbolic.einstein_tensor import einstein_tensor, field_equation
b = bg.symbolic_background(); H = b["H"]
print("a(t) =", b["a_t"], "   a(eta) =", b["a_eta"], "   calH =", b["calH"])
print("residuals of d eta/dt = 1/a and a(eta(t)) = exp(Ht):", bg.verify_time_relation())
"""),
code("""
rows = []
for name, key, coords in (("cosmic", "g_cosmic", b["coords_cosmic"]), ("conformal", "g_conformal", b["coords_conformal"])):
    G = Geometry(coords, b[key])
    res = G.constant_curvature_residual(H)
    rows.append((name, G.ricci_scalar, G.ricci_squared(), G.kretschmann(),
                 all(v == 0 for v in res.values()), field_equation(G, 3*H**2) == sp.zeros(4,4)))
for r in rows: print(r)
"""),
md("Sign convention test: with the Riemann sign flipped the Ricci scalar becomes $-12H^2$; the test suite treats that as a *detected* error."),
code("""
from desitter_perturbations.visualization import geometry_plots as gp, background_plots as bp
show(gp.plot_curvature_invariants())
show(bp.plot_background_evolution())
show(bp.plot_horizon_crossing())
"""),
md("""
Embedding check: flat slicing is the hyperboloid $-X_0^2+\\mathbf X^2+X_4^2=1/H^2$ with
$X_0=\\sinh(Ht)/H+\\tfrac H2 e^{Ht}r^2$, $X_4=\\cosh(Ht)/H-\\tfrac H2 e^{Ht}r^2$, $X_i=e^{Ht}x_i$. The induced metric is verified symbolically in
`validation.symbolic_checks.background_checks` (record `bg.embedding`).
"""),
])

NOTEBOOKS["02_metric_variation"] = ("Metric variation: inverse metric and determinant", [
md("""
# 02. Metric variation

Verify $g^{\\mu\\nu}=\\bar g^{\\mu\\nu}-\\epsilon h^{\\mu\\nu}+\\epsilon^2 h^\\mu{}_\\rho h^{\\rho\\nu}+O(\\epsilon^3)$ and
$\\sqrt{-g}=\\sqrt{-\\bar g}[1+\\tfrac\\epsilon2 h+\\epsilon^2(\\tfrac18 h^2-\\tfrac14 h_{\\mu\\nu}h^{\\mu\\nu})]$ with a **generic**
symmetric $h$ (ten independent symbols) and then show that the checks detect plausible mistakes.
"""),
code("""
from desitter_perturbations.symbolic import metric_expansion as me
gbar = me.generic_diagonal_background(); h = me.generic_perturbation()
print("inverse metric residual (should be zero matrix):")
print(me.check_inverse_metric(gbar, h, 2))
print("sqrt(-g) coefficient residuals:", me.check_sqrt_det(gbar, h))
"""),
code("""
print("wrong sign at 2nd order detected:", any(e != 0 for e in me.check_inverse_metric(gbar, h, 2, builder=me.inverse_metric_series_wrong_sign)))
d1, d2 = me.check_sqrt_det(gbar, h, formula=me.sqrt_det_series_formula_wrong_factor)
print("factor-of-two error in h_mn h^mn term detected:", d1 == 0 and d2 != 0)
print("residual of the wrong determinant formula:", d2)
"""),
md("Numerical order check: truncating the series one order earlier must degrade the error by exactly one power of $\\epsilon$."),
code("""
from desitter_perturbations.visualization import geometry_plots as gp
path, slopes = gp.plot_expansion_errors()
show(path); print({k: round(v, 3) for k, v in slopes.items()})
"""),
])

NOTEBOOKS["03_linearized_einstein_equations"] = ("Linearized Einstein equations and Bianchi identities", [
md("""
# 03. Linearized Einstein equations and Bianchi identities

Two independent derivations of $\\delta E_{\\mu\\nu}=\\delta G_{\\mu\\nu}+\\Lambda h_{\\mu\\nu}$ are compared for a generic polynomial perturbation:
(i) background-covariant formulas and (ii) the exact $O(\\epsilon)$ coefficient of $G+\\Lambda g$ computed from a truncated inverse metric.
The linearized Bianchi identity is then tested in three forms, including the *naive* divergence, which does not vanish.
"""),
code("""
from desitter_perturbations.background import symbolic_background
from desitter_perturbations.symbolic.curvature import Geometry
from desitter_perturbations.symbolic import einstein_tensor as et, bianchi_identities as bi
b = symbolic_background(); H = b["H"]; eta = b["eta"]; x, y, z = b["spatial"]
geom = Geometry(b["coords_conformal"], b["g_conformal"])
h = sp.Matrix([[eta*x, y, 0, z*eta],[y, x**2, eta*z, 0],[0, eta*z, y*z, x],[z*eta, 0, x, eta**2]])
cov = et.linearized_einstein_covariant(geom, h, 3*H**2)
e0, e1 = et.linearized_einstein_bruteforce(b["coords_conformal"], b["g_conformal"], h, 3*H**2)
print("background part of the brute-force pipeline (must vanish):", e0 == sp.zeros(4,4))
print("covariant minus component-level dE:", (cov["dE"] - e1).applyfunc(lambda e: sp.cancel(sp.together(e))) == sp.zeros(4,4))
"""),
code("""
res, naive, corr = bi.bianchi_covariant_residual(geom, cov["dG"], cov["dGamma"])
print("covariant-index Bianchi residual:", res)
print("mixed-index Bianchi residual   :", bi.bianchi_mixed_residual(geom, h, cov["dG"], cov["dGamma"]))
print("nabla^m dE_mn                  :", bi.bianchi_field_equation_residual(geom, cov["dE"]))
print()
print("naive nabla^m dG_mn (NOT zero; proportional to Lambda):")
for n in naive: print("   ", sp.factor(n))
"""),
md("""
The naive divergence is nonzero because $\\delta(\\nabla^\\mu G_{\\mu\\nu})$ contains $\\delta\\Gamma\\cdot\\bar G$ terms and
$\\bar G_{\\mu\\nu}=-\\Lambda\\bar g_{\\mu\\nu}\\neq0$. For $E=G+\\Lambda g$ the background value vanishes, so $\\bar\\nabla^\\mu\\delta E_{\\mu\\nu}=0$ holds exactly.
"""),
])

NOTEBOOKS["04_gauge_symmetry"] = ("Gauge symmetry", [
md("""
# 04. Gauge symmetry

Infinitesimal diffeomorphisms act as $h_{\\mu\\nu}\\to h_{\\mu\\nu}-\\bar\\nabla_\\mu\\xi_\\nu-\\bar\\nabla_\\nu\\xi_\\mu$.
We verify (i) $\\mathcal L_\\xi\\bar g=\\bar\\nabla\\xi+\\bar\\nabla\\xi$, (ii) that pure-gauge perturbations solve the linearized equation,
(iii) the scalar-vector transformation laws and residual gauge freedom, and (iv) a numerical comparison of one configuration in three gauges.
"""),
code("""
from desitter_perturbations.background import symbolic_background
from desitter_perturbations.symbolic.curvature import Geometry
from desitter_perturbations.symbolic import gauge_transformations as gt, einstein_tensor as et
b = symbolic_background(); H=b["H"]; eta=b["eta"]; x,y,z=b["spatial"]
geom = Geometry(b["coords_conformal"], b["g_conformal"])
xi = [eta**2*x*y, x*z**2 + eta, y*eta*z, z**3*x]
print("Lie derivative equals covariant symmetrised gradient:", gt.lie_equals_covariant_residual(xi, geom) == sp.zeros(4,4))
hpg = gt.pure_gauge_perturbation(xi, geom)
print("dE[h_pure_gauge] == 0:", et.linearized_einstein_covariant(geom, hpg, 3*H**2)["dE"] == sp.zeros(4,4))
"""),
code("""
laws = gt.svt_transformation_laws()
print("shift of A    :", laws["dA"]); print("shift of psi  :", laws["dpsi"])
print("shift of S_x  :", laws["dS_x"]); print("shift of B_z  :", laws["dB_z"])
print("residuals against closed-form laws:", gt.check_svt_laws())
"""),
md("""
Residual gauge freedom: Newtonian gauge ($B=E=0$) keeps $T=L=0$ for $k\\neq0$ (none). Synchronous gauge ($A=B=0$) keeps $T=C/a$ (one function of space).
Spatially flat gauge ($\\psi=E=0$) requires $\\mathcal H T=0$, hence none for $\\mathcal H\\neq0$.
"""),
code("""
from desitter_perturbations.visualization import gauge_plots as gaup
show(gaup.plot_gauge_comparison())
"""),
])

NOTEBOOKS["05_scalar_vector_tensor_decomposition"] = ("Scalar-vector-tensor decomposition", [
md("""
# 05. Scalar-vector-tensor decomposition

Algebraic identities are proven for a fully generic symmetric $h_{ij}(\\mathbf k)$; the numerical decomposition on a periodic box is then
tested for reconstruction error, transversality, resolution dependence, the $\\mathbf k=0$ mode and small-$k$ error amplification.
"""),
code("""
from desitter_perturbations.symbolic import tensor_decomposition as td
res = td.verify_svt_identities()
for k, v in td.all_zero(res).items(): print(f"{k:24s} {v}")
"""),
code("""
from desitter_perturbations.perturbations import svt_numeric as sv
for N in (8, 12, 16, 24):
    print(N, "round-trip error", sv.reconstruction_error(sv.random_symmetric_field(N, seed=N)))
g = sv.tt_part_field(16, seed=4)
print("spectral div  :", np.abs(sv.spectral_divergence(g)).max())
print("2nd-order FD  :", np.abs(sv.finite_difference_divergence(g)).max(), "(discretisation matters)")
"""),
code("""
from desitter_perturbations.visualization import gauge_plots as gaup
show(gaup.plot_svt_diagnostics())
"""),
md("""
Boundary conditions: the SVT split is a Helmholtz-type decomposition; on a periodic box it is unique except for the homogeneous mode,
which is assigned to $\\psi_0$ (trace) and a constant traceless $\\gamma_0$. On non-periodic domains uniqueness depends on falloff conditions and none of these routines apply.
"""),
])

NOTEBOOKS["06_quadratic_action"] = ("Second-order action", [
md("""
# 06. Second-order Einstein-Hilbert action

The tensor action is derived two ways and compared up to total derivatives (Euler-Lagrange operator of the difference):
route 1 expands $\\sqrt{-g}(R-2\\Lambda)$ exactly in $\\epsilon$; route 2 uses $L_2=-\\tfrac12\\sqrt{-\\bar g}\\,h^{\\mu\\nu}\\delta E_{\\mu\\nu}$.
Negative controls: a doubled normalisation and a wrong $\\Lambda$ are both detected.
"""),
code("""
from desitter_perturbations.symbolic import quadratic_action as qa
s = qa.setup()
r1 = qa.route1_direct_expansion(); r2 = qa.route2_from_linearized_equations()
print("L0 =", r1["L0"], "   L1 =", r1["L1"])
tgt = qa.target_tensor_lagrangian(s)
print("target (a^2/4)(gamma'^2 - (d gamma)^2) =", sp.simplify(tgt))
"""),
code("""
print("route 1 - target, E-L:", qa.total_derivative_check(sp.simplify(r1["L2"] - tgt), s))
print("route 2 - target, E-L:", qa.total_derivative_check(sp.simplify(r2["L2"] - tgt), s))
print("route 1 - route 2, E-L:", qa.total_derivative_check(sp.simplify(r1["L2"] - r2["L2"]), s))
print("doubled normalisation (must be non-zero):", qa.total_derivative_check(sp.simplify(r2["L2"] - 2*tgt), s))
"""),
code("""
fo = qa.first_order_variation_check()
print("first-order variation, Lambda = 3H^2, E-L expressions:", fo["correct"]["EL"])
print("first-order variation, Lambda = 2H^2, E-L expressions:", fo["wrong"]["EL"])
print("equation of motion from the action:", qa.tensor_equation_from_action(s)[0])
"""),
md("""
Parametrisation caveat: the code uses $g_{ij}=a^2(\\delta_{ij}+\\epsilon\\gamma_{ij})$, whereas many references use $a^2(e^\\gamma)_{ij}$. They differ at $O(\\gamma^2)$ by a term
proportional to the background field equations, so the quadratic actions agree up to boundary terms (see `docs/quadratic_action.md`).
"""),
])

NOTEBOOKS["07_tensor_mode_dynamics"] = ("Tensor-mode dynamics", [
md("""
# 07. Tensor-mode dynamics

Equation: $\\gamma_k''+2\\mathcal H\\gamma_k'+k^2\\gamma_k=0$ with exact Bunch-Davies solution $\\gamma_k=\\frac{iH}{M_{pl}k^{3/2}}(1+ik\\eta)e^{-ik\\eta}$.
We integrate in four formulations (conformal time, cosmic time, e-folds, canonical variable) and compare each against the exact solution and against each other.
"""),
code("""
from desitter_perturbations.dynamics import convergence as cv, residuals as rs
sols, errs = cv.formulation_comparison(k=1.0, rtol=1e-10)
for f, e in errs.items(): print(f"{f:10s} max relative error vs exact: {e:.2e}   nfev = {sols[f].nfev}")
print("cross-formulation max differences:", {k: f"{v:.2e}" for k, v in rs.cross_formulation(sols).items()})
"""),
code("""
from desitter_perturbations.visualization import mode_plots as mp, phase_space_plots as pp
show(mp.plot_mode_evolution()); show(mp.plot_formulation_comparison()); show(pp.plot_phase_space())
"""),
md("""
Physical reading: while $x=k/(aH)\\gg1$ the mode oscillates and its amplitude falls as $1/a$ (the physical wavelength is inside the Hubble radius);
for $x\\ll1$ the gradient term is negligible, $\\gamma_k\\to$ const, and the amplitude is frozen. Scalar-field mode equations are used only as a mathematical analogy.
"""),
])

NOTEBOOKS["08_horizon_crossing"] = ("Horizon crossing and numerical error", [
md("""
# 08. Horizon crossing

The crossover $x=1$ is where subhorizon and superhorizon approximations both fail. Check the exact amplitude there ($\\sqrt2$ in units of $H/(M_{pl}k^{3/2})$),
confirm there is no error spike, and use residual diagnostics that do not rely on knowing the exact solution.
"""),
code("""
from desitter_perturbations.spectra import mode_spectrum as msp
print("|gamma| at x = 1 in reference units:", msp.horizon_crossing_amplitude(), "  exact sqrt(2) =", np.sqrt(2))
"""),
code("""
from desitter_perturbations.dynamics import convergence as cv
traj = cv.error_trajectories(rtols=(1e-4, 1e-7, 1e-10))
for rtol, (x, e) in traj.items():
    near = e[(x > 0.3) & (x < 3)].max(); late = e[x < 0.05].max()
    print(f"rtol {rtol:.0e}: max error near x=1 {near:.2e}, late-time plateau {late:.2e}")
"""),
code("""
from desitter_perturbations.visualization import mode_plots as mp, convergence_plots as cp
show(cp.plot_tolerance_and_trajectory()); show(mp.plot_residual_diagnostics())
"""),
md("""
The residual diagnostic is limited by interpolation of the stored output (it converges as the 4th power of output density), not by the integrator.
"""),
])

NOTEBOOKS["09_gauge_invariant_diagnostics"] = ("Gauge-invariant diagnostics and constraints", [
md("""
# 09. Gauge-invariant diagnostics and constraints

In a general gauge the scalar-sector linearized equations depend on $(A,B,\\psi,E)$ only through $\\Phi_B$ and $\\Psi_B$; the vector sector only through $V=S-F'$.
The vacuum constraints then remove the scalar and vector sectors for $k\\neq0$, leaving the two tensor polarizations.
"""),
code("""
from desitter_perturbations.perturbations import scalar_sector, vector_sector, tensor_sector, constraints
print("scalar residuals :", scalar_sector.general_gauge_dE()["residuals"])
print("vector residuals :", vector_sector.vector_dE()["residuals"])
t = tensor_sector.tensor_dE()["residuals"]; print("tensor residuals :", {k: t[k] for k in ("xx","yy","xy")}, t["others"])
"""),
code("""
c = constraints.scalar_constraint_solution()
print("traceless ij  :", c["e_traceless"], "  => Phi = Psi")
print("0i            :", c["e0z"])
print("00            :", c["e00"])
print("00 after constraints:", c["e00_after_constraints"], "  => Psi = 0 for k != 0")
"""),
md("Counting: 10 metric components, 4 constraints, 4 gauge parameters give 2 propagating degrees of freedom."),
code("""
r = tensor_sector.covariant_tensor_equation()
print("transverse:", r["divergence"], " traceless:", r["trace"], " dE + (box h - 2H^2 h)/2 = 0:", r["equation"] == sp.zeros(4,4))
"""),
])

NOTEBOOKS["10_canonical_normalization"] = ("Canonical normalization", [
md("""
# 10. Canonical normalization

With $e^\\lambda_{ij}e^{\\lambda'}_{ij}=2\\delta^{\\lambda\\lambda'}$ the tensor action becomes $S=\\sum_\\lambda\\frac12\\int z_T^2[\\gamma_\\lambda'^2-(\\partial\\gamma_\\lambda)^2]$.
The prefactor $z_T$ is *derived* here from the symbolic Lagrangian of notebook 06.
"""),
code("""
from desitter_perturbations.quantization import canonical_variables as cv, mode_normalization as mn
z2, d = cv.derive_z_T_squared(); print("z_T^2 =", z2)
print("z''/z =", cv.effective_potential())
el, zpp = cv.verify_canonical_action(); print("E-L of L[gamma=v/z] - L_v (zero => total derivative):", el)
eq, e, v, zc = cv.mode_equation_from_action(); print("canonical equation of motion:", eq, "= 0")
"""),
code("""
C2, W = mn.normalization_from_wronskian()
print("Wronskian condition fixes |C|^2 =", C2)
print("exact mode residual:", mn.verify_mode_solves_equation(), "   positive-frequency limit:", mn.positive_frequency_limit())
"""),
code("""
from desitter_perturbations.visualization import phase_space_plots as pp
show(pp.plot_phase_space())
"""),
])

NOTEBOOKS["11_tensor_power_spectrum"] = ("Tensor power spectrum", [
md("""
# 11. Tensor power spectrum

Definition: $\\langle\\gamma_{ij}\\gamma_{ij}\\rangle=\\int d\\ln k\\,\\Delta_t^2$ with $\\Delta_t^2=\\frac{n_{pol}N_e k^3}{2\\pi^2}|\\gamma_k|^2$, $n_{pol}=2$, $N_e=2$.
For the exact Bunch-Davies mode this gives $\\Delta_t^2=\\frac{2H^2}{\\pi^2M_{pl}^2}(1+k^2\\eta^2)$, which is checked against numerically evolved modes.
"""),
code("""
from desitter_perturbations.spectra import polarization as pol, tensor_power as tp, parameter_sweeps as ps
k = np.array([0.3, -1.1, 0.8])
print("normalisation error  :", pol.check_normalization(k))
print("completeness error   :", pol.check_completeness(k))
print("transverse/traceless :", pol.check_transverse_traceless(k))
r = ps.sweep_wavenumber(n=12, k_min=0.2, k_max=100, eta_end=-1e-2, x_start=300)
print("max relative error of numerical spectrum:", r["rel_err"].max())
"""),
code("""
from desitter_perturbations.visualization import spectrum_plots as sp_
show(sp_.plot_tensor_spectrum()); show(sp_.plot_late_time())
"""),
md("The amplitude $2H^2/(\\pi^2M_{pl}^2)$ is a benchmark for this implementation under the stated conventions; it is not a statement about the observed universe."),
])

NOTEBOOKS["12_parameter_sweeps"] = ("Parameter sweeps", [
md("""
# 12. Parameter sweeps

Sweeps over $H$ and $M_{pl}$ (scaling $\\Delta_t^2\\propto H^2/M_{pl}^2$), over Bogoliubov-rotated initial states, and over the start time of the leading-WKB initial condition.
The last is an *exploratory observation with an exact explanation*: leading-WKB initial data at finite $x_{start}$ is a Bogoliubov state with $|\\beta|=1/(2x_{start}^2)$,
and the numerical power ratio equals the exact prediction to $\\sim10^{-10}$.
"""),
code("""
from desitter_perturbations.spectra import parameter_sweeps as ps
r = ps.sweep_H_Mpl(); print("max |numeric/predicted - 1|:", np.abs(r["numeric"]/r["predicted"] - 1).max())
b = ps.sweep_bogoliubov(); print("max |numeric - predicted| (Bogoliubov):", np.abs(b["numeric"] - b["predicted"]).max())
rows = ps.sweep_x_start()
for xs, ratio, pred, beta in zip(rows["x_start"], rows["ratio"], rows["predicted"], rows["abs_beta"]):
    print(f"x_start {xs:6.0f}  P/P_BD {ratio:.10f}  predicted {pred:.10f}  |beta| {beta:.3e}  1/(2x^2) {1/(2*xs**2):.3e}")
"""),
code("""
from desitter_perturbations.visualization import spectrum_plots as sp_
show(sp_.plot_bogoliubov_and_start()); show(sp_.plot_scaling_H_M())
"""),
code("""
import csv, os
path = os.path.join("..", "data", "generated", "sweep_x_start.csv")
print(open(path).read() if os.path.exists(path) else "run scripts/run_parameter_sweeps.py first")
"""),
])

NOTEBOOKS["13_numerical_convergence"] = ("Numerical convergence", [
md("""
# 13. Numerical convergence

Tolerance study, solver comparison, fixed-step RK4 order of convergence (expected 4), and a 30-digit mpmath cross-check of the double-precision solver.
"""),
code("""
from desitter_perturbations.dynamics import convergence as cv
for r in cv.tolerance_study(): print(f"rtol {r['rtol']:.0e}  max err {r['max_err']:.2e}  nfev {r['nfev']}")
h, e, p = cv.rk4_order_study(); print("RK4 observed order:", round(p, 3))
for r in cv.method_comparison(): print(f"{r['method']:7s} err {r['max_err']:.2e}  nfev {r['nfev']:6d}  ok {r['success']}")
"""),
code("""
from desitter_perturbations.dynamics.numerical_solvers import solve_mode, solve_mode_mpmath, SolverConfig
xe, g_mp = solve_mode_mpmath(1.0, 10.0, 0.1, dps=30)
s = solve_mode(1.0, 10.0, 0.1, "conformal", SolverConfig(rtol=1e-12, atol=1e-15, n_out=3))
print("mpmath:", g_mp); print("DOP853:", s.gamma[-1]); print("relative difference:", abs(s.gamma[-1]-g_mp)/abs(g_mp))
"""),
code("""
from desitter_perturbations.visualization import convergence_plots as cp
show(cp.plot_methods_and_order())
"""),
])

NOTEBOOKS["14_full_reproducibility_study"] = ("Full reproducibility study", [
md("""
# 14. Full reproducibility study

Runs the complete verification suite (symbolic identities, analytic benchmarks, numerical checks), prints the summary and shows how each check is recorded.
This takes a few minutes because the symbolic checks are exact computations.
"""),
code("""
from desitter_perturbations.validation import validation_report as vr
recs, meta = vr.run_all(fast=False, symbolic=True)
print(meta)
print(len(recs), "records;", sum(r.passed for r in recs), "passed;", sum(r.kind == "negative_control" for r in recs), "negative controls")
"""),
code("""
import collections
by = collections.defaultdict(list)
for r in recs: by[r.category].append(r)
for c, rs in by.items(): print(f"{c:22s} {len(rs):3d} records   all passed: {all(r.passed for r in rs)}")
"""),
code("""
worst = sorted([r for r in recs if r.kind in ("numerical", "benchmark")], key=lambda r: -r.residual / max(r.tolerance, 1e-300))[:8]
print("checks closest to their tolerance (residual / tolerance):")
for r in worst: print(f"  {r.id:36s} {r.residual:.2e} / {r.tolerance:.1e}")
"""),
md("""
What the report establishes and what it does not: it shows that the implemented equations are mutually consistent, reproduce independent analytic results and are numerically converged;
it does not show that the free-field de Sitter model is a complete description of gravity (see `docs/limitations.md`).
"""),
])


def build(name, title, cells):
    nb = nbf.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.metadata["title"] = title
    nb.cells.append(nbf.v4.new_code_cell(SETUP))
    for kind, text in cells:
        nb.cells.append(nbf.v4.new_markdown_cell(text) if kind == "md" else nbf.v4.new_code_cell(text))
    return nb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-execute", action="store_true")
    ap.add_argument("--only", nargs="*", default=None)
    a = ap.parse_args()
    NB.mkdir(exist_ok=True)
    for name, (title, cells) in NOTEBOOKS.items():
        if a.only and name[:2] not in a.only:
            continue
        nb = build(name, title, cells)
        if not a.no_execute:
            NotebookClient(nb, timeout=1500, kernel_name="python3", resources={"metadata": {"path": str(NB)}}).execute()
        nbf.write(nb, NB / f"{name}.ipynb")
        print("wrote", NB / f"{name}.ipynb", flush=True)


if __name__ == "__main__":
    main()
