"""Numerical integration of tensor modes in four formulations.

The user chooses the formulation, solver, tolerances, start/end of the integration in
the horizon variable x = k/(aH) = -k eta (x decreases with time), the initial data
and the output density (logarithmic in x).  Results are returned in a common
container so that formulations can be compared point by point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
from scipy.integrate import solve_ivp

from . import analytic_solutions as an
from . import initial_conditions as ic
from . import tensor_modes as tm


@dataclass
class SolverConfig:
    method: str = "DOP853"
    rtol: float = 1e-10
    atol: float = 1e-14          # absolute tolerance in units of the reference amplitude H/(M k^1.5)
    n_out: int = 400             # number of (logarithmically spaced in x) output points
    max_step: float = np.inf
    H: float = 1.0
    M_pl: float = 1.0


@dataclass
class ModeSolution:
    formulation: str
    k: float
    x: np.ndarray                 # horizon variable at output points (decreasing)
    eta: np.ndarray
    t: np.ndarray                 # cosmic time with a(t=0)=1
    gamma: np.ndarray             # complex gamma_k
    dgamma_deta: np.ndarray       # complex d gamma_k / d eta
    nfev: int = 0
    success: bool = True
    message: str = ""
    config: Optional[SolverConfig] = None
    ic_name: str = ""

    @property
    def a(self):
        return -1.0 / (self.config.H * self.eta)


def _initial_state(ic_name, k, eta0, cfg: SolverConfig, bogoliubov=None):
    if ic_name == "bd_exact":
        return ic.bd_exact(k, eta0, cfg.H, cfg.M_pl)
    if ic_name == "bd_asymptotic":
        return ic.bd_asymptotic(k, eta0, cfg.H, cfg.M_pl)
    if ic_name == "bogoliubov":
        alpha, beta = bogoliubov
        return ic.bogoliubov(k, eta0, alpha, beta, cfg.H, cfg.M_pl)
    raise ValueError(f"unknown initial condition {ic_name!r}")


def solve_mode(k: float, x_start: float = 100.0, x_end: float = 1e-3,
               formulation: str = "conformal", config: Optional[SolverConfig] = None,
               ic_name: str = "bd_exact", bogoliubov=None) -> ModeSolution:
    """Integrate one tensor mode from x_start (early, subhorizon) to x_end (late).

    Parameters
    ----------
    k : comoving wavenumber (same units as H, with a(t=0) = 1 in the cosmic formulations).
    formulation : 'conformal' | 'cosmic' | 'efold' | 'canonical'.
    """
    cfg = config or SolverConfig()
    if formulation not in tm.FORMULATIONS:
        raise ValueError(f"formulation must be one of {tm.FORMULATIONS}")
    if not (x_start > x_end > 0):
        raise ValueError("require x_start > x_end > 0")
    H, M = cfg.H, cfg.M_pl
    eta0 = -x_start / k
    g0, dg0 = _initial_state(ic_name, k, eta0, cfg, bogoliubov)
    x_out = np.geomspace(x_start, x_end, cfg.n_out)
    eta_out = -x_out / k
    gref = an.gamma_ref(k, H, M)
    atol = cfg.atol * gref
    opts = dict(method=cfg.method, rtol=cfg.rtol, atol=atol, max_step=cfg.max_step)

    if formulation == "conformal":
        y0 = tm.pack(g0, dg0)
        sol = solve_ivp(tm.rhs_conformal, (eta_out[0], eta_out[-1]), y0, t_eval=eta_out, args=(k,), **opts)
        g, dg = tm.unpack(sol.y)
    elif formulation == "cosmic":
        a0 = -1.0 / (H * eta0)
        t0 = np.log(a0) / H
        t_out = -np.log(-H * eta_out) / H
        # derivative wrt t:  gdot = g'/a
        y0 = tm.pack(g0, dg0 / a0)
        sol = solve_ivp(tm.rhs_cosmic, (t_out[0], t_out[-1]), y0, t_eval=t_out, args=(k, H), **opts)
        g, gdot = tm.unpack(sol.y)
        dg = gdot * np.exp(H * t_out)
    elif formulation == "efold":
        a0 = -1.0 / (H * eta0)
        N0 = np.log(a0)
        N_out = np.log(-1.0 / (H * eta_out))
        y0 = tm.pack(g0, dg0 / (a0 * H))   # g_N = gdot/H = g'/(a H)
        sol = solve_ivp(tm.rhs_efold, (N_out[0], N_out[-1]), y0, t_eval=N_out, args=(k, H), **opts)
        g, gN = tm.unpack(sol.y)
        dg = gN * H * np.exp(N_out)
    else:  # canonical
        z0 = an.z_T(eta0, H, M)
        calH0 = -1.0 / eta0
        v0 = z0 * g0
        dv0 = z0 * (calH0 * g0 + dg0)
        y0 = tm.pack(v0, dv0)
        sol = solve_ivp(tm.rhs_canonical, (eta_out[0], eta_out[-1]), y0, t_eval=eta_out, args=(k,), **opts)
        v, dv = tm.unpack(sol.y)
        z = an.z_T(eta_out, H, M)
        g = v / z
        dg = (dv - (-1.0 / eta_out) * v) / z

    t_cos = -np.log(-H * eta_out) / H
    return ModeSolution(formulation=formulation, k=float(k), x=x_out, eta=eta_out, t=t_cos,
                        gamma=g, dgamma_deta=dg, nfev=int(sol.nfev), success=bool(sol.success),
                        message=str(sol.message), config=cfg, ic_name=ic_name)


# ------------------------------------------------------------------ fixed-step RK4

def solve_rk4_efold(k: float, x_start: float, x_end: float, n_steps: int,
                    H: float = 1.0, M_pl: float = 1.0, ic_name: str = "bd_exact"):
    """Classical RK4 with a fixed step in N = ln a.  Used for order-of-convergence studies.

    Returns (x_final, gamma_final, dgamma_deta_final).
    """
    eta0 = -x_start / k
    g0, dg0 = _initial_state(ic_name, k, eta0, SolverConfig(H=H, M_pl=M_pl))
    a0 = -1.0 / (H * eta0)
    N0 = np.log(a0)
    N1 = N0 + np.log(x_start / x_end)
    h = (N1 - N0) / n_steps
    y = np.array([g0, dg0 / (a0 * H)], dtype=complex)

    def f(N, y):
        return np.array([y[1], -3.0 * y[1] - (k * k * np.exp(-2.0 * N) / H**2) * y[0]])

    N = N0
    for _ in range(n_steps):
        k1 = f(N, y)
        k2 = f(N + h / 2, y + h / 2 * k1)
        k3 = f(N + h / 2, y + h / 2 * k2)
        k4 = f(N + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        N += h
    a1 = np.exp(N)
    eta1 = -1.0 / (H * a1)
    return -k * eta1, y[0], y[1] * H * a1


# --------------------------------------------------------- arbitrary precision

def solve_mode_mpmath(k: float, x_start: float = 20.0, x_end: float = 1e-2, dps: int = 30,
                      H: float = 1.0, M_pl: float = 1.0):
    """High-precision Taylor-series integration in conformal time with mpmath.odefun.

    Only practical for modest x_start (the cost grows with the number of oscillations);
    used as an independent cross-check of the double-precision solvers.
    Returns (x_end, gamma_end) as Python complex numbers.
    """
    import mpmath as mp

    mp.mp.dps = dps
    kk = mp.mpf(k)
    eta0 = -mp.mpf(x_start) / kk
    gref = mp.mpf(H) / (mp.mpf(M_pl) * kk ** mp.mpf(1.5))
    g0 = 1j * gref * (1 + 1j * kk * eta0) * mp.e ** (-1j * kk * eta0)
    dg0 = 1j * gref * kk**2 * eta0 * mp.e ** (-1j * kk * eta0)

    def F(eta, y):
        return [y[1], (2 / eta) * y[1] - kk**2 * y[0]]

    sol = mp.odefun(F, eta0, [g0, dg0])
    eta1 = -mp.mpf(x_end) / kk
    val = sol(eta1)[0]
    return float(x_end), complex(val)
